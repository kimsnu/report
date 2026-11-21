#!/usr/bin/env python3
"""Find closed reciprocal-best-hit cycles across an ordered list of species.

Python 3.8+. DIAMOND is the default search engine for genome-scale data.
The small-data Smith-Waterman engine uses only the Python standard library.

    python3 ring_rbh.py --type cds --outdir result A=A.fa B=B.fa C=C.fa
    python3 ring_rbh.py --manifest species.tsv --threads 32 --outdir result

A manifest has a header and three tab-separated columns (species, type, fasta).
Type is aa or cds. Relative paths are resolved against the manifest directory.
List the first species only once. Each FASTA ID must be unique within its species.

CDS are translated in frame 1 with genetic code 1. Stop at the first definite
stop and ignore 1-2 trailing bases for the protein search. A degenerate codon
encoding a single amino acid is resolved, and other ambiguous codons become X.
For protein inputs, use the part before the first '*'. Search sequences are
uppercase, with U and O represented by X. Original CDS are retained for output.

For each adjacent species pair, rank targets independently in both directions.
Rank the best local alignment of each target by raw score among eligible hits.
Tied best targets are discarded by default. With --ties lexicographic, the
smallest target ID is selected.
Keep a group only if every link is reciprocal and the last link returns to the
same first-species record. This tests the edges of the specified ring only.

Outputs include ortholog_sets.tsv (original FASTA IDs), rbh_pairs.tsv (edge
evidence), run.json (settings and counts), and sets/OG*.aa.fasta. Each set also has an
OG*.cds.fasta with original nucleotides when every input species supplied CDS.
FASTA headers start with species|original_ID and retain the source description.
For protein inputs the output FASTA retains the original amino acid sequence.
"""

import argparse
import csv
import gzip
import itertools
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import NamedTuple, Optional


def log(message):
    print(message, file=sys.stderr, flush=True)


class Source(NamedTuple):
    name: str
    kind: str
    path: Path


class Record(NamedTuple):
    seq_id: str
    header: str
    protein: str
    original: str


def get_sources(args):
    if bool(args.manifest) == bool(args.inputs):
        raise ValueError("supply either ordered SPECIES=FASTA arguments or --manifest FILE")
    sources = []
    if args.manifest:
        with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
            table = csv.reader(handle, delimiter="\t")
            if next(table, None) != ["species", "type", "fasta"]:
                raise ValueError("manifest header must be: species<TAB>type<TAB>fasta")
            for number, row in enumerate(table, 2):
                if not row or all(not field.strip() for field in row):
                    continue
                if len(row) != 3:
                    raise ValueError(f"manifest line {number}: expected three columns")
                name, kind, path = (value.strip() for value in row)
                sources.append(Source(name, kind, args.manifest.parent / path))
    else:
        for item in args.inputs:
            name, sep, filename = item.partition("=")
            if not sep or not filename:
                raise ValueError(f"expected SPECIES=FASTA, received {item!r}")
            sources.append(Source(name, args.type, Path(filename)))
    names = set()
    for name, kind, path in sources:
        if not name or any(ch.isspace() or ch in "|=" for ch in name):
            raise ValueError(f"invalid species label {name!r}: use a label without whitespace, '|' or '='")
        if name in names:
            raise ValueError(f"duplicate species {name!r}: the ring closes automatically")
        if kind not in {"aa", "cds", "auto"}:
            raise ValueError(f"{name}: type must be aa or cds (or auto for tagged inputs)")
        if not path.is_file():
            raise ValueError(f"{name}: input file does not exist: {path}")
        names.add(name)
    if len(sources) < 2:
        raise ValueError("at least two different species are required")
    return sources


def fasta_records(path):
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as stream:
        header, chunks = None, []
        for number, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
                if not header:
                    raise ValueError(f"{path}:{number}: empty FASTA header")
            elif header is None:
                raise ValueError(f"{path}:{number}: expected a FASTA header")
            else:
                chunks.extend(line.split())
        if header is None:
            raise ValueError(f"{path}: empty FASTA")
        yield header, "".join(chunks)


def detect_type(header, path):
    tokens = header.split(None, 2)
    tag = tokens[1].split(":")[0].lower() if len(tokens) > 1 else ""
    if tag == "cds":
        return "cds"
    if tag in {"pep", "peptide", "protein", "aa"}:
        return "aa"
    tags = set(path.name.lower().split("."))
    aa, cds = bool(tags & {"faa", "pep", "peptide", "protein", "aa"}), bool(tags & {"cds", "ffn"})
    if aa != cds:
        return "aa" if aa else "cds"
    raise ValueError(f"{path}: specify --type aa or --type cds, or put the type in a manifest")


IUPAC = dict(zip("ACGT", "ACGT"))
IUPAC.update(R="AG", Y="CT", S="GC", W="AT", K="GT", M="AC", B="CGT",
             D="AGT", H="ACT", V="ACG", N="ACGT")
CODE = dict(zip(("".join(c) for c in itertools.product("TCAG", repeat=3)),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
DNA_CHARS = set(IUPAC) | {"U"}
AA_CHARS = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*")
SEARCH_AA = str.maketrans({"U": "X", "O": "X"})


@lru_cache(maxsize=4096)
def ambiguous_codon(codon):
    choices = {CODE["".join(b)] for b in itertools.product(*(IUPAC[c] for c in codon))}
    return next(iter(choices)) if len(choices) == 1 else "X"


def translate(seq):
    bases = seq.replace("U", "T")
    peptide = []
    for i in range(0, len(bases) - 2, 3):
        codon = bases[i:i + 3]
        aa = CODE.get(codon) or ambiguous_codon(codon)
        if aa == "*":
            break
        peptide.append(aa)
    return "".join(peptide)


def prepare_source(src):
    kind = src.kind
    records, ids = [], set()
    for header, original in fasta_records(src.path):
        seq_id = header.split(None, 1)[0]
        if seq_id in ids:
            raise ValueError(f"{src.name}: duplicate FASTA ID {seq_id!r}")
        ids.add(seq_id)
        if kind == "auto":
            kind = detect_type(header, src.path)
        seq = original.upper()
        invalid = set(seq) - (DNA_CHARS if kind == "cds" else AA_CHARS)
        if invalid or not seq:
            raise ValueError(f"{src.name}/{seq_id}: invalid or empty {kind} sequence (unexpected symbols {sorted(invalid)})")
        protein = translate(seq) if kind == "cds" else seq.split("*", 1)[0].translate(SEARCH_AA)
        if protein:
            records.append(Record(seq_id, header, protein, original))
    if not records:
        raise ValueError(f"{src.name}: no nonempty proteins after translation or stop handling")
    return kind, records, len(ids) - len(records)


class Hit(NamedTuple):
    score: int
    identity: float
    qcov: float
    scov: float
    evalue: Optional[float] = None


@dataclass
class Best:
    target: int
    hit: Hit
    tied: bool = False


def consider(best, qi, si, hit, targets, args):
    if (hit.score < args.min_score or hit.identity < args.min_id or
            hit.qcov < args.min_qcov or hit.scov < args.min_scov):
        return
    previous = best.get(qi)
    if previous is None or hit.score > previous.hit.score:
        best[qi] = Best(si, hit)
    elif hit.score == previous.hit.score and si != previous.target:
        previous.tied = True
        if targets[si].seq_id < targets[previous.target].seq_id:
            previous.target, previous.hit = si, hit


def reciprocal(forward, reverse, ties):
    edges = {}
    for qi, f in forward.items():
        r = reverse.get(f.target)
        if r is not None and r.target == qi and (ties == "lexicographic" or not (f.tied or r.tied)):
            edges[qi] = (f.target, f, r)
    return edges


def write_fasta(stream, header, seq):
    stream.write(f">{header}\n")
    for i in range(0, len(seq), 60):
        stream.write(seq[i:i + 60] + "\n")


def make_databases(species, work, exe, threads):
    databases = []
    for i, records in enumerate(species):
        query = work / f"species{i}.faa"
        db = work / f"species{i}"
        with query.open("w", encoding="ascii", newline="\n") as stream:
            for j, rec in enumerate(records):
                # Internal IDs preserve arbitrary source accessions, pipes and version suffixes.
                write_fasta(stream, f"r{j}", rec.protein)
        subprocess.run([exe, "makedb", "--in", str(query), "--db", str(db),
                        "--threads", str(threads), "--quiet"], check=True)
        databases.append((query, db))
    return databases


def diamond_best(queries, targets, query_file, database, exe, work, args):
    fields = "qseqid sseqid score nident length qstart qend sstart send evalue".split()
    command = [exe, "blastp", "--query", str(query_file), "--db", str(database),
               "--outfmt", "6", *fields, "--max-target-seqs", "0", "--max-hsps", "1",
               "--evalue", str(args.evalue), "--threads", str(args.threads),
               "--sensitive", "--tmpdir", str(work), "--quiet"]
    best = {}
    with subprocess.Popen(command, stdout=subprocess.PIPE, encoding="ascii") as process:
        try:
            for line in process.stdout:
                q, s, score, identity_count, length, q1, q2, s1, s2, evalue = line.split()
                qi, si = int(q[1:]), int(s[1:])
                if (q != f"r{qi}" or s != f"r{si}" or not 0 <= qi < len(queries)
                        or not 0 <= si < len(targets)):
                    raise ValueError(f"unexpected DIAMOND sequence IDs: {q} {s}")
                hit = Hit(int(score), 100 * int(identity_count) / int(length),
                          100 * (abs(int(q2) - int(q1)) + 1) / len(queries[qi].protein),
                          100 * (abs(int(s2) - int(s1)) + 1) / len(targets[si].protein), float(evalue))
                consider(best, qi, si, hit, targets, args)
        except BaseException:
            process.kill()
            raise
        status = process.wait()
        if status:
            raise subprocess.CalledProcessError(status, command)
    return best


# NCBI BLOSUM62 lower triangle. X has score -1 against every residue.
# https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62
def blosum62():
    alphabet = b"ARNDCQEGHILKMFPSTWYVBJZ"
    triangle = """
    4
    -1 5
    -2 0 6
    -2 -2 1 6
    0 -3 -3 -3 9
    -1 1 0 0 -3 5
    -1 0 0 2 -4 2 5
    0 -2 0 -1 -3 -2 -2 6
    -2 0 1 -1 -3 0 0 -2 8
    -1 -3 -3 -3 -1 -3 -3 -4 -3 4
    -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4
    -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5
    -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5
    -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6
    -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7
    1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4
    0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5
    -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11
    -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7
    0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4
    -2 -1 4 4 -3 0 1 -1 0 -3 -4 0 -3 -3 -2 0 -1 -4 -3 -3 4
    -1 -2 -3 -3 -1 -2 -3 -4 -3 3 3 -3 2 0 -3 -2 -1 -2 -1 2 -3 3
    -1 0 0 1 -3 4 4 -2 0 -3 -3 1 -1 -3 -1 0 -1 -2 -2 -2 0 -3 4
    """.strip().splitlines()
    matrix = [[-1] * 91 for _ in range(91)]
    for i, row in enumerate(triangle):
        for j, value in enumerate(row.split()):
            x, y = alphabet[i], alphabet[j]
            matrix[x][y] = matrix[y][x] = int(value)
    return matrix


MATRIX = blosum62()


def smith_waterman(query, target):
    """Local affine alignment. A gap of k residues costs 11 + (k - 1)."""
    a, b = query.encode("ascii"), target.encode("ascii")
    width = len(b) + 1
    trace = bytearray((len(a) + 1) * width)
    prev, vertical = [0] * width, [-1000000] * width
    best, end_i, end_j = 0, 0, 0
    for i, x in enumerate(a, 1):
        cur, horizontal = [0] * width, -1000000
        scores, row_start = MATRIX[x], i * width
        for j, y in enumerate(b, 1):
            up_open, up_extend = prev[j] - 11, vertical[j] - 1
            left_open, left_extend = cur[j - 1] - 11, horizontal - 1
            flags = 0
            if up_extend > up_open:
                vertical[j], flags = up_extend, 4
            else:
                vertical[j] = up_open
            if left_extend > left_open:
                horizontal, flags = left_extend, flags | 8
            else:
                horizontal = left_open
            score, move = prev[j - 1] + scores[y], 1
            if vertical[j] > score:
                score, move = vertical[j], 2
            if horizontal > score:
                score, move = horizontal, 3
            if score <= 0:
                score, move = 0, 0
            cur[j], trace[row_start + j] = score, flags | move
            if score > best:
                best, end_i, end_j = score, i, j
        prev = cur
    if best == 0:
        return None
    i, j = end_i, end_j
    state = identical = columns = 0
    while True:
        flags = trace[i * width + j]
        move = state or (flags & 3)
        if move == 0:
            break
        columns += 1
        if move == 1:
            i, j, state = i - 1, j - 1, 0
            identical += a[i] == b[j]
        elif move == 2:
            i, state = i - 1, (2 if flags & 4 else 0)
        else:
            j, state = j - 1, (3 if flags & 8 else 0)
    return Hit(best, 100 * identical / columns, 100 * (end_i - i) / len(a),
               100 * (end_j - j) / len(b))


def python_best(left, right, args):
    forward, reverse = {}, {}
    for i, q in enumerate(left):
        for j, s in enumerate(right):
            hit = smith_waterman(q.protein, s.protein)
            if hit is not None:
                consider(forward, i, j, hit, right, args)
                flipped = Hit(hit.score, hit.identity, hit.scov, hit.qcov)
                consider(reverse, j, i, flipped, left, args)
    return forward, reverse


def closed_cycles(links, starting_records):
    for first in sorted(links[0], key=lambda i: starting_records[i].seq_id):
        members, current = [first], first
        for link in links:
            current = link.get(current)
            if current is None:
                break
            members.append(current)
        else:
            if current == first:
                yield members[:-1]


def evidence(best):
    h = best.hit
    return [h.score, "NA" if h.evalue is None else f"{h.evalue:g}", f"{h.identity:.3f}",
            f"{h.qcov:.3f}", f"{h.scov:.3f}", int(best.tied)]


def search_ring(sources, species, dbs, exe, work, output, args):
    links, counts = [], []
    n = len(sources)
    pairs = 1 if n == 2 else n
    with (output / "rbh_pairs.tsv").open("w", encoding="utf-8", newline="") as file:
        table = csv.writer(file, delimiter="\t", lineterminator="\n")
        metrics = ["raw_score", "evalue", "identity_pct", "qcov_pct", "scov_pct", "tied"]
        table.writerow(["from_species", "from_id", "to_species", "to_id"] +
                       [d + "_" + m for d in ("forward", "reverse") for m in metrics])
        for i in range(pairs):
            j = (i + 1) % n
            a, b = sources[i].name, sources[j].name
            log(f"[{i + 1}/{pairs}] {a} <-> {b}")
            if args.engine == "python":
                forward, reverse = python_best(species[i], species[j], args)
            else:
                forward = diamond_best(species[i], species[j], dbs[i][0], dbs[j][1], exe, work, args)
                reverse = diamond_best(species[j], species[i], dbs[j][0], dbs[i][1], exe, work, args)
            edges = reciprocal(forward, reverse, args.ties)
            mapping = {}
            for q, (s, f, r) in edges.items():
                mapping[q] = s
                table.writerow([a, species[i][q].seq_id, b, species[j][s].seq_id] + evidence(f) + evidence(r))
            links.append(mapping)
            if n == 2:
                links.append({s: q for q, s in mapping.items()})
            counts.append(dict(source=a, target=b, reciprocal_pairs=len(edges),
                               forward_top_ties=sum(h.tied for h in forward.values()),
                               reverse_top_ties=sum(h.tied for h in reverse.values())))
            log(f"  {len(edges):,} reciprocal pairs")
    return links, counts


def write_sets(output, sources, species, kinds, links):
    folder = output / "sets"
    folder.mkdir()
    all_cds = all(kind == "cds" for kind in kinds)
    count = 0
    with (output / "ortholog_sets.tsv").open("w", encoding="utf-8", newline="") as file:
        table = csv.writer(file, delimiter="\t", lineterminator="\n")
        table.writerow(["set_id"] + [s.name for s in sources])
        for count, members in enumerate(closed_cycles(links, species[0]), 1):
            group = f"OG{count:06d}"
            records = [seqs[i] for seqs, i in zip(species, members)]
            table.writerow([group] + [r.seq_id for r in records])
            for ext in (("aa", "cds") if all_cds else ("aa",)):
                with (folder / f"{group}.{ext}.fasta").open("w", encoding="utf-8", newline="\n") as fa:
                    for src, kind, rec in zip(sources, kinds, records):
                        sequence = rec.original if ext == "cds" or kind == "aa" else rec.protein
                        write_fasta(fa, f"{src.name}|{rec.header}", sequence)
    return count


def run(args):
    sources = get_sources(args)
    out = args.outdir.absolute()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"output directory must be new or empty: {out}")
    exe = version = None
    if args.engine == "diamond":
        exe = shutil.which(args.diamond)
        if not exe:
            raise ValueError("DIAMOND not found. Install DIAMOND, supply --diamond PATH, "
                             "or choose --engine python for small datasets.")
        version = subprocess.check_output([exe, "version"], text=True).strip()
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".ring-rbh-", dir=out.parent) as tmp:
        work = Path(tmp)
        output = work / "result"
        output.mkdir()
        names = [s.name for s in sources]
        log("Ring: " + " -> ".join(names + names[:1]))
        species, kinds, inputs = [], [], []
        for src in sources:
            kind, records, skipped = prepare_source(src)
            species.append(records)
            kinds.append(kind)
            inputs.append(dict(species=src.name, fasta=str(src.path.resolve()), type=kind,
                               searchable=len(records), empty_translation=skipped))
            log(f"{src.name} ({kind}): {len(records):,} searchable, {skipped} empty translations")
        dbs = make_databases(species, work, exe, args.threads) if exe else None
        links, counts = search_ring(sources, species, dbs, exe, work, output, args)
        groups = write_sets(output, sources, species, kinds, links)
        meta = dict(engine=args.engine, diamond_version=version,
                    threads=args.threads if exe else 1, genetic_code=1,
                    diamond_sensitivity="sensitive" if exe else None, min_raw_score=args.min_score,
                    min_identity_pct=args.min_id, min_qcov_pct=args.min_qcov, min_scov_pct=args.min_scov,
                    max_evalue=args.evalue, ties=args.ties, inputs=inputs, adjacent_pairs=counts,
                    closed_sets=groups)
        (output / "run.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        # Publish a completed directory. Failed runs cannot leave stale set FASTAs.
        if out.exists():
            out.rmdir()  # Only an empty directory can be removed.
        output.rename(out)
    log(f"Closed sets: {groups:,} | {out}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        epilog="Defaults are DIAMOND sensitive, E-value 1e-5, at least 50% coverage of each sequence, "
               "and exclusion of tied best targets. The Python engine does not estimate E-values.")
    p.add_argument("inputs", nargs="*", metavar="SPECIES=FASTA", help="species in ring order")
    p.add_argument("--manifest", type=Path, help="ordered TSV with species, type, fasta columns")
    p.add_argument("--type", choices=["auto", "aa", "cds"], default="auto", help="sequence type for positional files")
    p.add_argument("--outdir", type=Path, default=Path("ring_rbh_results"), help="new or empty output directory (default %(default)s)")
    p.add_argument("--engine", choices=["diamond", "python"], default="diamond",
                   help="DIAMOND sensitive blastp for proteomes, or exhaustive Smith-Waterman for small datasets")
    p.add_argument("--diamond", default="diamond", help="DIAMOND executable name or path")
    p.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1), help="threads per DIAMOND process (default %(default)s)")
    p.add_argument("--evalue", type=float, help="maximum E-value for DIAMOND (default 1e-5)")
    p.add_argument("--min-score", type=int, default=0, help="minimum raw alignment score")
    p.add_argument("--min-id", type=float, default=0, help="minimum identical alignment columns (percent)")
    p.add_argument("--min-qcov", type=float, default=50, help="minimum query coverage (percent)")
    p.add_argument("--min-scov", type=float, default=50, help="minimum subject coverage (percent)")
    p.add_argument("--ties", choices=["drop", "lexicographic"], default="drop", help="policy for equal best target scores")
    args = p.parse_args()
    if any(not 0 <= v <= 100 for v in (args.min_id, args.min_qcov, args.min_scov)):
        p.error("identity and coverage percentages must be in [0, 100]")
    if args.min_score < 0 or args.threads < 1:
        p.error("--min-score must be nonnegative and --threads must be positive")
    if args.evalue is not None:
        if args.engine != "diamond" or not math.isfinite(args.evalue) or args.evalue <= 0:
            p.error("--evalue requires the DIAMOND engine and a finite positive value")
    elif args.engine == "diamond":
        args.evalue = 1e-5
    try:
        run(args)
    except (ValueError, OSError, EOFError, csv.Error, subprocess.SubprocessError) as error:
        p.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()

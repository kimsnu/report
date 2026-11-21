#!/usr/bin/env python3
"""Keep the longest translated isoform of each gene in each input FASTA.

Python 3.8+, standard library only. Inputs may be plain text or gzip.
Each input produces a .longest.fasta sequence file and a .longest.tsv selection table.

    python3 longest_translation.py --type cds --outdir longest_cds "*.cds.all.fa"
    python3 longest_translation.py --type aa --outdir longest_aa "*.pep.all.fa.gz"

CDS must be spliced, oriented 5' to 3', and in frame from the first base.
Use the standard genetic code. Count complete codons before the first in-frame
stop. A partially ambiguous codon counts as one residue unless it always
encodes a stop. For protein input count residues before the first '*'.
The chosen sequence is written in its original alphabet and case.
"""

import argparse
import csv
import glob
import gzip
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple


GENE_PATTERN = r"(?:^|[\s\[])gene[:=]([^\s;\]]+)"
TRANSCRIPT = re.compile(r"(?:^|[\s\[])transcript[:=]([^\s;\]]+)")
INVALID_CDS = re.compile(r"[^ACGTURYSWKMBDHVN]")
INVALID_AA = re.compile(r"[^A-Z*]")
# TAR expands to TAA/TAG and TRA expands to TAA/TGA.
STOP_CODONS = re.compile(r"TAA|TAG|TGA|TAR|TRA")
FASTA_SUFFIXES = {".fa", ".fas", ".fasta", ".fna", ".ffn", ".faa", ".fsa"}


def fasta_records(path):
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        header = None
        parts = []
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_no}: empty FASTA header")
                parts = []
            else:
                if header is None:
                    raise ValueError(f"{path}:{line_no}: expected a FASTA header starting with '>'")
                parts.extend(line.split())
        if header is None:
            raise ValueError(f"{path}: no FASTA records")
        yield header, "".join(parts)


def detect_type(header, path):
    fields = header.split(None, 2)
    marker = fields[1].split(":", 1)[0].lower() if len(fields) > 1 else ""
    if marker == "cds":
        return "cds"
    if marker in {"pep", "peptide", "protein", "aa"}:
        return "aa"
    labels = set(path.name.lower().split("."))
    is_cds = bool(labels & {"cds", "ffn"})
    is_aa = bool(labels & {"aa", "pep", "peptide", "protein", "faa"})
    if is_cds != is_aa:
        return "cds" if is_cds else "aa"
    raise ValueError(f"{path}: sequence type unclear. Supply --type cds or --type aa.")


def translated_length(sequence, kind):
    if not sequence:
        raise ValueError("empty sequence")
    seq = sequence.upper()
    invalid = (INVALID_CDS if kind == "cds" else INVALID_AA).search(seq)
    if invalid:
        raise ValueError(f"invalid symbol {invalid[0]!r} at position {invalid.start() + 1} in {kind} sequence")
    if kind == "aa":
        stop = seq.find("*")
        return len(seq) if stop == -1 else stop
    seq = seq.replace("U", "T")
    for stop in STOP_CODONS.finditer(seq):
        if stop.start() % 3 == 0:
            return stop.start() // 3
    return len(seq) // 3


@dataclass
class Isoform:
    header: str
    sequence: str
    aa_length: int
    # Prefer a shared transcript ID for ties when Ensembl supplies one.
    tie_key: Tuple[str, str]  # (transcript ID or FASTA ID, FASTA ID)


@dataclass
class Selection:
    best: Isoform
    n_isoforms: int = 1
    n_tied_longest: int = 1

    def consider(self, candidate):
        self.n_isoforms += 1
        if candidate.aa_length > self.best.aa_length:
            self.best = candidate
            self.n_tied_longest = 1
        elif candidate.aa_length == self.best.aa_length:
            self.n_tied_longest += 1
            if candidate.tie_key < self.best.tie_key:
                self.best = candidate


def select_isoforms(path, kind, gene_pattern):
    genes = {}
    for header, sequence in fasta_records(path):
        if kind == "auto":
            kind = detect_type(header, path)
        seq_id = header.split(None, 1)[0]
        match = gene_pattern.search(header)
        if not match or not match.group(1):
            raise ValueError(f"{path}: {seq_id}: cannot find gene ID. Check the header or use --gene-regex.")
        gene_id = match.group(1)
        transcript = TRANSCRIPT.search(header)
        tie_id = transcript.group(1) if transcript else seq_id
        try:
            length = translated_length(sequence, kind)
        except ValueError as error:
            raise ValueError(f"{path}: {seq_id}: {error}") from error
        candidate = Isoform(header, sequence, length, (tie_id, seq_id))
        if gene_id in genes:
            genes[gene_id].consider(candidate)
        else:
            genes[gene_id] = Selection(candidate)
    return kind, genes


def write_outputs(genes, fasta_path, tsv_path):
    # Stage both outputs before replacing any previous run for this file.
    with tempfile.TemporaryDirectory(dir=fasta_path.parent) as directory:
        staged_fasta = Path(directory) / fasta_path.name
        staged_tsv = Path(directory) / tsv_path.name
        with staged_fasta.open("w", encoding="utf-8", newline="\n") as fasta, \
                staged_tsv.open("w", encoding="utf-8", newline="") as table:
            rows = csv.writer(table, delimiter="\t", lineterminator="\n")
            rows.writerow(["gene_id", "sequence_id", "aa_length", "n_isoforms", "n_tied_longest"])
            for gene_id, selection in genes.items():
                best = selection.best
                fasta.write(f">{best.header}\n")
                for start in range(0, len(best.sequence), 60):
                    fasta.write(best.sequence[start:start + 60] + "\n")
                rows.writerow([gene_id, best.tie_key[1], best.aa_length,
                               selection.n_isoforms, selection.n_tied_longest])
        staged_fasta.replace(fasta_path)
        staged_tsv.replace(tsv_path)


def plan_outputs(patterns, outdir):
    inputs = []
    for pattern in patterns:
        path = Path(pattern)
        matches = [path] if path.is_file() else sorted(Path(p) for p in glob.glob(pattern) if Path(p).is_file())
        if not matches:
            raise ValueError(f"no input file matches {pattern!r}")
        inputs.extend(matches)
    jobs = []
    used = {path.resolve() for path in inputs}
    for path in inputs:
        name = Path(path.stem) if path.suffix.lower() == ".gz" else Path(path.name)
        base = name.stem if name.suffix.lower() in FASTA_SUFFIXES else name.name
        fasta = outdir / f"{base}.longest.fasta"
        tsv = outdir / f"{base}.longest.tsv"
        for target in (fasta, tsv):
            if target.resolve() in used:
                raise ValueError(f"file name collision at {target}. Use distinct input basenames or output directories.")
            if target.exists() and not target.is_file():
                raise ValueError(f"output path is not a file: {target}")
            used.add(target.resolve())
        jobs.append((path, fasta, tsv))
    return jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("inputs", nargs="+", metavar="FASTA", help="input FASTA files or glob patterns")
    parser.add_argument("--type", choices=("auto", "cds", "aa"), default="auto",
                        help="one sequence type per file; auto uses the FASTA header or filename (default)")
    parser.add_argument("--outdir", type=Path, default=Path("longest"),
                        help="directory for *.longest.fasta sequences and *.longest.tsv selection tables (default: longest)")
    parser.add_argument("--gene-regex", default=GENE_PATTERN, metavar="REGEX",
                        help="regex with one capture group for the gene ID, applied to the header without '>' "
                             "(default: gene:ID or gene=ID)")
    args = parser.parse_args()
    try:
        gene_pattern = re.compile(args.gene_regex)
        if gene_pattern.groups != 1:
            raise ValueError("--gene-regex must contain exactly one capture group for the gene ID")
        jobs = plan_outputs(args.inputs, args.outdir)
        args.outdir.mkdir(parents=True, exist_ok=True)
        for path, fasta_path, tsv_path in jobs:
            kind, genes = select_isoforms(path, args.type, gene_pattern)
            write_outputs(genes, fasta_path, tsv_path)
            total = sum(selection.n_isoforms for selection in genes.values())
            print(f"{path.name} ({kind}): {total:,} records -> {len(genes):,} genes | {fasta_path} | {tsv_path}",
                  file=sys.stderr, flush=True)
            del genes  # Release the previous species before loading the next.
    except (OSError, EOFError, ValueError, re.error) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()

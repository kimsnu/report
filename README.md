# reporti
#REPresentative ORThologs Identifier

The main concept is designating a representative gene of in-paralogs (co-orthologs) which stands for the most conserved genes that are assigned based on the phylogenetic tree between all member of orthologs including the co-orthologs in the other species (Figure a). If the in-paralogs are ancestral genes, the software designates a representative gene among the ancestral in-paralogs and finally the descendant genes are assigned to the representative genes (Figure b).

<img width="942" height="778" alt="image" src="https://github.com/user-attachments/assets/437415d0-4825-4710-911a-d2853f1c55f8" />



# To make the longest translation
**CDS fasta**
python3 longest_translation.py --type cds --outdir longest_cds "*.fa"

**AA fasta**
python3 longest_translation.py --type aa --outdir longest_aa "*.fa"

# To make representative orthologs
sudo apt install diamond-aligner

python3 ring_rbh.py --type cds --threads 100 --outdir rbh_cds \\
  at=Arabidopsis_thaliana.TAIR10.cds.all.longest.fasta \\
  gm=Glycine_max.Glycine_max_v2.1.cds.all.longest.fasta \\
  mt=Medicago_truncatula.MtrunA17r5.0_ANR.cds.all.longest.fasta \\
  pv=Phaseolus_vulgaris.PhaVulg1_0.cds.all.longest.fasta \\
  vf=Vicia_faba.Hedin2_genome_v1.cds.all.longest.fasta

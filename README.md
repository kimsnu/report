# reporti
# REPresentative ORThologs Identifier

The main concept is designating a representative gene of in-paralogs (co-orthologs) which stands for the most conserved genes that are assigned based on the phylogenetic tree between all member of orthologs including the co-orthologs in the other species. If the in-paralogs are ancestral genes, the software designates a representative gene among the ancestral in-paralogs and finally the descendant genes are assigned to the representative genes.

<img width="463" height="389" alt="image" src="https://github.com/user-attachments/assets/55f5b840-2af1-4dd8-be96-50a6d643fa86" />

Figure. Many-to-many orthologs by speciation and gene duplication events. X and Y are duplication events, Z is a speciation node. α3A and α5B are representative orthologs.

# To make the longest translation
**CDS fasta**

python3 longest_translation.py --type cds --outdir longest_cds "*.fasta"

or

python3 longest_translation.py --type cds --outdir longest_cds "*.fasta.gz"


**AA fasta**

python3 longest_translation.py --type aa --outdir longest_aa "*.fasta"

or

python3 longest_translation.py --type aa --outdir longest_aa "*.fasta.gz"


# To make representative orthologs
sudo apt install diamond-aligner


python3 ring_rbh.py --type cds --threads 100 --outdir rbh_cds \\

  at=Arabidopsis_thaliana.TAIR10.cds.all.longest.fasta \\
  
  gm=Glycine_max.Glycine_max_v2.1.cds.all.longest.fasta \\
  
  mt=Medicago_truncatula.MtrunA17r5.0_ANR.cds.all.longest.fasta \\
  
  pv=Phaseolus_vulgaris.PhaVulg1_0.cds.all.longest.fasta \\
  
  vf=Vicia_faba.Hedin2_genome_v1.cds.all.longest.fasta

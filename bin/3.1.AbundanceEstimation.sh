##This script performs index generation of the human transcriptome and gene pseudomapping using Samlmon
##Run salmon index to generate the index. 
salmon index -t ../transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz hsa_index
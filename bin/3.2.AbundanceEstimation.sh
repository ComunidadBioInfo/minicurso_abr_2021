##This script performs pseudoalignment and abundance estimation of transcripts using salmon
mkdir -p ../salmon_quants
for i in {1..6};  do
r1="s${i}_R1.fastq.gz"
r2="s${i}_R2.fastq.gz"
printf "\n"
echo "Processing sample s$i"
##Use salmon to quantificate transcripts
~/salmon-latest_linux_x86_64/bin/salmon quant \
                                         -i ../transcriptome/hsa_index -l A \ #Path to the index
                                         -1 ../data/paired/$r1 \ #Path to r1 reads
                                         -2 ../data/paired/$r2 \ #Path to r2 reads
                                         -p 8 --validateMappings -o ../salmon_quants/s${i}_quant
done
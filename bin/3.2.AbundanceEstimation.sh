##This script performs pseudoalignment and abundance estimation of transcripts using salmon
mkdir -p ../salmon_quants
for i in {1..12};  do
r1="s${i}_R1.fastq.gz"
r2="s${i}_R2.fastq.gz"
printf "\n"
echo "Processing sample s$i"
##Use salmon to quantificate transcripts
~/salmon-latest_linux_x86_64/bin/salmon quant -i ../transcriptome/hsa_index -l A -1 ../data/$r1 -2 ../data/$r2 -p 8 --validateMappings -o ../salmon_quants/s${i}_quant
done

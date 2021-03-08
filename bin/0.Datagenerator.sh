mkdir -p split

## Split fastq files in chunks of 2M reads
for i in {1..12}; do
for j in {1..2}; do
zcat t${i}_R${j}.fastq.gz | split -l 8000000 - t${i}_R${j}.fastq
done
done

mv *.fastqa* split/
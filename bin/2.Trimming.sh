##This script performs trimming of reads at the 5' end using trimmomatic
#Create a directory to move all the unpaired files
mkdir -p ../data/unpaired
mkdir -p ../data/paired
#Set the input and output paths
in_path="../data"
output_paired_path="../data/paired"
output_un_path="../data/unpaired"

#Create arrays with the name of the input (forward and reverse) files
for i in {1..6}; do
printf "\n" #Insert a line break
echo "s$i" #Print the name of each sample
r1="s${i}_R1.fastq.gz"
r2="s${i}_R2.fastq.gz"
#Create arrays with the name of the output (paired and unpaired files)
r1p="s${i}_R1.paired.fastq.gz"
r2p="s${i}_R2.paired.fastq.gz"
r1u="s${i}_R1.unpaired.fastq.gz"
r2u="s${i}_R2.unpaired.fastq.gz"
#Run trimmomatic with HEADCROP
trimmomatic PE -threads 8 $in_path/$r1 $in_path/$r2 $output_paired_path/$r1p $output_un_path/$r1u $output_paired_path/$r2p $output_un_path/$r2u HEADCROP:2
done 
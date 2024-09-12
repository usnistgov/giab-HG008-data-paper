#Author: Lesley Chapman Hannah
#Date: August 2021
#Description:
#Script adds sets of annotations from ANNOVAR database. Update protocol and operation flags to add/remove annotations
#Input: vcf file
#Output: vcf with ANNOVAR annotations

#!/bin/bash
#SBATCH --time=5-0:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=200g
#SBATCH --job-name run_annovar
#SBATCH -o %j.out
#SBATCH -e %j.err

#Load modules
module load annovar
#Generate annovar vcf
table_annovar.pl $1 /fdb/annovar/current/hg38 --thread 10 --buildver hg38 --outfile $2 --remove --protocol refGene,spliceai_filtered,spidex,cosmic92_coding,cosmic92_noncoding,gnomad211_exome,gnomad211_genome,esp6500siv2_all,1000g2015aug_all,clinvar_20221231,clinvar_20220320 --operation g,f,f,f,f,f,f,f,f,f,f --vcfinput

#Add clinvar annotations
#wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20231230.vcf.gz .


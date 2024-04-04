#!/bin/bash -l
#SBATCH -t 14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marlo072@umn.edu
#SBATCH --job-name ANNOVAR
#SBATCH -o %j.ANNOVAR.out
#SBATCH -e %j.ANNOVAR.err
#SBATCH -p msismall,msilarge

conda activate snakeprac

perl table_annovar.pl ../../datafiles/SmallSample.vcf.gz equCab/ --vcfinput --outfile outputs/Chr1ANNOVAROutputTest2 --buildver equCab3 --protocol ensGene --operation g --codingarg '--tolerate' --nastring '.' --nopolish

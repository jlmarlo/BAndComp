#!/bin/bash -l
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marlo072@umn.edu
#SBATCH --job-name Annovar
#SBATCH -o %j.Annovar.out
#SBATCH -e %j.Annovar.err
#SBATCH -p msismall,msilarge


#MAW=Make Annovar Work
## This script contains the code used to generate the files necessary for ANNOVAR to run
#perl annotate_variation.pl --downdb --buildver equCab3 ensGene equCab/

##either use this code to download a reference from UCSC or just have your own reference nearby
#perl annotate_variation.pl --downdb --buildver equCab3 seq equCab/equCab3_seq/

perl retrieve_seq_from_fasta.pl equCab/equCab3_ensGene.txt -seqfile /panfs/jay/groups/6/durwa004/marlo072/Documents/ReferenceFiles/goldenPath.Ec_build-3.0_wMSY.fa -format ensGene \
	-outfile equCab/equCab3_ensGeneMrna.fa

##name of outfile must end with ensGeneMrna.fa as that is what the next command will look for
##Command with make a fasta file that contains the sequences found at the locations indicated by the gene definition file

## These are the files you need for set up. Code for actual annotation will be found in the pipeline 

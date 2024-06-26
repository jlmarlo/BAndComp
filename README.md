# BAndComp (Better Annotation Comparison)

This workflow takes a VCF and annotates all variants present with three annotator programs: Ensembl Variant Effect Predictor (VEP), SnpEff, and ANNOVAR

Annotations are compared with a custom python script and report about annotation agreement is generated in the output directory. 

The following instructions will help you download this repository, create an environment with necessary programs and annotate a blank VCF. 

## Setup
**Create environment with dependecies**
Assuming you already have conda or mamba installed

```
#create environment
mamba create -n BandCompEnv -c bioconda -c conda-forge snakemake=7.19.1 bcftools=1.17 gatk4=4.2.6.1 snpeff=5.2

```
**Download repository**

```
git clone https://github.com/jlmarlo/BAndComp.git
```
You will be prompted to enter your github username and a password
**The password is not your github password. It is an access token**
To create an access token:
1. Click on your profile avatar in the top right of the screen
2. Click on Settings
3. On the left side at the bottom of the list of options click on Developer settings
4. Click Personal Access Tokens
5. Click Tokens (classic)
6. Click Generate new token at the top
7. Click Generate new token (classic)
8. In the Note box type whatever you'd like
9. Select the box next to repo
10. Scroll to bottom and click Generate token
11. Copy the access token. It's a series of random letters and numbers. And save it somewhere.

You can copy and paste this token any time you are asked for a github password.

**Download singularity image for VEP usage**

```
cd BAndComp
wget https://s3.msi.umn.edu/wags/wags.sif
mkdir Jobs
mkdir datafiles
mv wags.sif datafiles/wags.sif
```

**Move datafiles into proper locations**

```
mkdir Jobs
mkdir datafiles
cp /home/durwa004/shared/AnnotationFiles/* datafiles/
cp /home/durwa004/shared/joint_call.goldenPath.20230726.vcf.gz* datafiles/
```

**Edit slurm script**

```
vim MarchBang.slurm
```
Change the '--mail-user' option to include your email

## Running program

**Check program integrity**
```
mamba activate BandCompEnv
snakemake -s BandComp.smk -np --configfile config.yaml

```
This should create an output that lists all steps that need to be completed. There should be no warning text at the bottom of this text it should finish with a table of numbers and step names. 

```
sbatch MarchBang.slurm
```

That should be running :)

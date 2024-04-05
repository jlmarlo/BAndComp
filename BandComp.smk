###This will be the snakemake file that contains the rules to create a VCF with all three annotations
##BAndComp Better Annotation Comparison

import os

#singularity: config['sif']

chromosomes = []
for x in range(1,32):
	x=str(x)
	chrom = 'chr'+x
	chromosomes.append(chrom)
chromosomes.append('chrX')
chromosomes.append('MSY')
chromosomes.append('chrM')


##rule all
## final combined vcf with all three annotations
rule all:
	input: 'outputs/final/PPID3Annotation.vcf.gz'

##rule split
##read in large VCF and split into chromosomes

rule split:
	input:
		vcf = 'datafiles/PPID_20230726.vcf.gz',
		tbi = 'datafiles/PPID_20230726.vcf.gz.tbi'
	output:
		smallVCF= 'outputs/splitVCFs/{chrom}/{chrom}_noAnno.vcf.gz',
		smalltbi= 'outputs/splitVCFs/{chrom}/{chrom}_noAnno.vcf.gz.tbi'
	resources:
		mem_mb=50000,
		time=420
	shell:
		'''
			bcftools view {input.vcf} -r {wildcards.chrom} \
				-Oz -o {output.smallVCF} 
			gatk IndexFeatureFile -I {output.smallVCF}
		'''

##rule vep
##take in a chromosome VCF and annotation with VEP

rule veping:
	input:
		vcf= 'outputs/splitVCFs/{chrom}/{chrom}_noAnno.vcf.gz',
		tbi= 'outputs/splitVCFs/{chrom}/{chrom}_noAnno.vcf.gz.tbi'
	output:
		vepTbi= 'outputs/Veped/{chrom}/{chrom}_Veped.vcf.gz.tbi',
		vepVCF = 'outputs/Veped/{chrom}/{chrom}_Veped.vcf.gz'
	params:
		ref_fasta= 'datafiles/goldenPath.Ec_build-3.0_wMSY.fa',
		gtf = 'datafiles/Equus_caballus.EquCab3.0.104_sorted.gtf.gz',
		unzipped = 'outputs/Veped/{chrom}/{chrom}_Veped.vcf'
	threads: 6
	singularity : config['sif']
	resources:
		time=420,
		mem_mb=60000
	shell:
		'''
			set -e
			source activate ensembl-vep

			vep --help

			vep -i {input.vcf} -o {output.vepVCF} \
				--fasta {params.ref_fasta} \
				--fork {threads} \
				--force_overwrite \
				--vcf \
				--gtf {params.gtf} \
				--dir_plugins /opt/wags/src/VepPlugins \
				--dont_skip \
				--protein \
				--variant_class \
				--biotype \
				--compress_output bgzip \
				--stats_text 
			#bgzip --threads {threads} -c {params.unzipped} > {output.vepVCF}
			tabix -p vcf {output.vepVCF}
		'''

##rule snpEff
##take in VCF annotated with VEP and annotate with snpEff

rule snpEff:
	input:
		vepvcf = 'outputs/Veped/{chrom}/{chrom}_Veped.vcf.gz',
		veptbi = 'outputs/Veped/{chrom}/{chrom}_Veped.vcf.gz.tbi'
	output:
		snpeffvcftbi = 'outputs/snpeff/{chrom}/{chrom}_SnpEff_Veped.vcf.gz.tbi',
		snpeffgz = 'outputs/snpeff/{chrom}/{chrom}_SnpEff_Veped.vcf.gz',
		stats = 'outputs/snpeff/summaries/{chrom}.stats.csv',
		html = 'outputs/snpeff/html/{chrom}.summary.html'
	params:
		splicerange = '2',
		updownrange = '5000',
		database = 'EquCab3.0.105',
		unzipped = 'outputs/snpeff/{chrom}/{chrom}_SnpEff_Veped.vcf'
	threads: 6
	resources:
		time=420,
		mem_mb=80000
	shell:
		'''
			snpEff -Xmx14g \
				{params.database} \
				-ss {params.splicerange} \
				-ud {params.updownrange} \
				{input.vepvcf} \
				> {params.unzipped} \
				-csvStats {output.stats} \
				-stats {output.html}
			bgzip --threads {threads} -c {params.unzipped} > {output.snpeffgz}
			gatk IndexFeatureFile -I {output.snpeffgz}
				
		'''

##rule annovar
##take in VCF annotated with VEP and snpEff and annotate with ANNOVAR
rule prep_annovar:
	output:
		ensGene = 'outputs/annovar_database/equCab3_ensGene.txt',
		mRNAFasta = 'outputs/annovar_database/equCab3_ensGeneMrna.fa'
	params:
		program = 'programs/annovar/',
		database = 'equCab3',
		outdir = 'outputs/annovar_database/',
		fasta = 'datafiles/goldenPath.Ec_build-3.0_wMSY.fa'
	resources:
		time=420,
		mem_mb=200000
	shell:
		'''
			perl {params.program}annotate_variation.pl --downdb \
				--buildver {params.database} \
				ensGene \
				{params.outdir}
			perl {params.program}retrieve_seq_from_fasta.pl {output.ensGene} \
				-seqfile {params.fasta} \
				-format ensGene \
				-outfile {output.mRNAFasta}
		'''


rule annovar:
	input:
		snpeffvcf = 'outputs/snpeff/{chrom}/{chrom}_SnpEff_Veped.vcf.gz',
		snpeffvcftbi = 'outputs/snpeff/{chrom}/{chrom}_SnpEff_Veped.vcf.gz.tbi',
		mRNAFasta = 'outputs/annovar_database/equCab3_ensGeneMrna.fa'
	output:
		annovarVCF = 'outputs/annovarVCFs/{chrom}/{chrom}_Annovar_SnpEff_Veped.vcf.gz',
		annovarVCFtbi = 'outputs/annovarVCFs/{chrom}/{chrom}_Annovar_SnpEff_Veped.vcf.gz.tbi',
		avInput = 'outputs/annovarVCFs/{chrom}/{chrom}_Annovar_SnpEff_Veped.avinput',
	params:
		program = 'programs/annovar/',
		datadirectory = 'outputs/annovar_database/',
		version = 'equCab3',
		outName = 'outputs/annovarVCFs/{chrom}/{chrom}_Annovar_SnpEff_Veped',
		outdir = 'outputs/annovarVCFs/{chrom}/',
		unzipped = 'outputs/annovarVCFs/{chrom}/{chrom}_Annovar_SnpEff_Veped.equCab3_multianno.vcf'
	resources:
		time=420,
		mem_mb=80000
	threads: 6
	shell:
		'''
			#mkdir {params.outdir}
			perl {params.program}table_annovar.pl {input.snpeffvcf} \
				{params.datadirectory} \
				--vcfinput \
				--outfile {params.outName} \
				--buildver {params.version} \
				--protocol ensGene \
				--operation g \
				--codingarg '--tolerate' \
				--nastring '.' \
				--nopolish
			bgzip --threads {threads} -c {params.unzipped} > {output.annovarVCF}
			gatk IndexFeatureFile -I {output.annovarVCF}
		'''

##rule combine
##take in all vcfs that have been annotated with all three annotators
##output vcf found in all rule
def gather_vcfs(wildcards):
	return expand('outputs/annovarVCFs/{chrom}/{chrom}_Annovar_SnpEff_Veped.vcf.gz',
		chrom=chromosomes)

rule combine:
	input: 
		gather_vcfs
	output:
		finalvcf = 'outputs/final/PPID3Annotation.vcf.gz',
		tbi = 'outputs/final/PPID3Annotation.vcf.gz.tbi'
	params:
		vcfs = lambda wildcards, input: " --input ".join(map(str,input))
	threads: 12
	resources:
		time =420,
		mem_mb=200000
	shell:
		'''
			gatk GatherVcfsCloud \
				--ignore-safety-checks \
				--gather-type BLOCK \
				--input {params.vcfs} \
				--output {output.finalvcf}
			gatk --java-options "-Xmx18g -Xms6g" \
				Index Feature File -I {output.finalvcf}
		'''

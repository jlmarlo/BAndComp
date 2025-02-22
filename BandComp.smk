###This will be the snakemake file that contains the rules to create a VCF with all three annotations
##BAndComp Better Annotation Comparison

##creates list of relevant chromosomes. Can be altered to include Y or mitochondrial if desired
chromosomes = []
for x in range(1,32):
	x=str(x)
	chrom = 'chr'+x
	chromosomes.append(chrom)
chromosomes.append('chrX')
#chromosomes.append('MSY')
#chromosomes.append('chrM')
VCF = config['VCF']
outputVCF = config['final']
##rule all
#
## final combined vcf with all three annotations
#And final annotation table
rule all:
	input: 
		expand('outputs/final/{outputVCF}', outputVCF=outputVCF),
		'outputs/final/AnnotationTable.tsv',
		'outputs/final/FinalReport.txt'



##rule split
##read in large VCF and split into chromosomes
###parallelizes workflow

rule split:
	input:
		vcf = expand('datafiles/{VCF}',VCF=VCF),
		tbi = expand('datafiles/{VCF}.tbi',VCF=VCF)
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

##decomposes multiallelic sites by splitting each allele onto a separate line of the VCF
rule decompose:
	input:
		vcf = 'outputs/splitVCFs/{chrom}/{chrom}_noAnno.vcf.gz',
		tbi = 'outputs/splitVCFs/{chrom}/{chrom}_noAnno.vcf.gz.tbi'
	output:
		decomp = 'outputs/decomposed/{chrom}/{chrom}_noAnno_decomp.vcf.gz',
		decomptbi = 'outputs/decomposed/{chrom}/{chrom}_noAnno_decomp.vcf.gz.tbi'
	resources:
		mem_mb = 50000,
		time = 420
	shell:
		'''
			bcftools norm {input.vcf} -m - -Oz -o {output.decomp}
			gatk IndexFeatureFile -I {output.decomp}
		'''


##rule vep
##take in a chromosome VCF and annotation with VEP
#Uses VEP version from Jonah's singularity image because I couldn't get VEP to work on my own
#Might be a sticking point

rule veping:
	input:
		vcf= 'outputs/decomposed/{chrom}/{chrom}_noAnno_decomp.vcf.gz',
		tbi= 'outputs/decomposed/{chrom}/{chrom}_noAnno_decomp.vcf.gz.tbi'
	output:
		vepTbi= 'outputs/Veped/{chrom}/{chrom}_Veped.vcf.gz.tbi',
		vepVCF = 'outputs/Veped/{chrom}/{chrom}_Veped.vcf.gz'
	params:
		ref_fasta= config['ref'],
		gtf = config['gtf'],
		unzipped = 'outputs/Veped/{chrom}/{chrom}_Veped.vcf',
		updown = config['updown']
	threads: 6
	singularity : config['sif']
	resources:
		time=420,
		mem_mb=60000
	shell:
		'''
			set -e
			source activate ensembl-vep

			vep -i {input.vcf} -o {output.vepVCF} \
				--fasta {params.ref_fasta} \
				--fork {threads} \
				--force_overwrite \
				--vcf \
				--gtf {params.gtf} \
				--dir_plugins /opt/wags/src/VepPlugins \
				--dont_skip \
				--protein \
				--distance {params.updown} \
				--variant_class \
				--biotype \
				--format vcf \
				--compress_output bgzip \
				--stats_text 
			tabix -p vcf {output.vepVCF}
		'''
                


##rule snpEff
##take in VCF annotated with VEP and annotate with snpEff
##This rule currently seems to be making files in the home directory of the program and not in
#it's individual output folder
rule snpEff:
	input:
		vepvcf = 'outputs/Veped/{chrom}/{chrom}_Veped.vcf.gz',
		veptbi = 'outputs/Veped/{chrom}/{chrom}_Veped.vcf.gz.tbi'
	output:
		snpeffvcftbi = 'outputs/snpeff/{chrom}/{chrom}_SnpEff_Veped.vcf.gz.tbi',
		snpeffgz = 'outputs/snpeff/{chrom}/{chrom}_SnpEff_Veped.vcf.gz',
		stats = 'outputs/snpeff/{chrom}/{chrom}.stats.csv',
		html = 'outputs/snpeff/{chrom}/{chrom}.summary.html'
	params:
		splicerange = '2',
		updownrange = config['updown'],
		database = config['snpEffdatabase'],
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


#ANNOVAR requires specific database set ups and files to be processed a specific way. This rule does this outside
#of the annovar annotation rule because it only needs to be run a single time. 
#currently set up only for equCab3 and the goldenPath reference. In the future 
#will need to be modified to potentially allow for other references. 


rule prep_annovar:
	input: 
		gtf = config['gtf']
	output:
		PreensGene = 'outputs/annovar_database/EquCab3_ensGene01.txt',
		ensGene = 'outputs/annovar_database/equCab3_ensGene02.txt',
		ensGeneFin = 'outputs/annovar_database/equCab3_ensGene.txt',
		mRNAFasta = 'outputs/annovar_database/equCab3_ensGeneMrna.fa'
	params:
		program = 'programs/annovar/',
		database = 'equCab3',
		outdir = 'outputs/annovar_database/',
		fasta = config['ref']
	resources:
		time=420,
		mem_mb=200000
	shell:
		'''
			gtfToGenePred -genePredExt {input.gtf} {output.PreensGene}
			nl -nln {output.PreensGene} > {output.ensGene}
			awk 'BEGIN{{OFS="\t"}} {{$3="chr"$3; print}}' {output.ensGene} > {output.ensGeneFin}
			perl {params.program}retrieve_seq_from_fasta.pl {output.ensGeneFin} \
				-seqfile {params.fasta} \
				-format ensGene \
				-outfile {output.mRNAFasta}
		'''

##This is the rule were annovar actually does the annotation on each chromosome VCF
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
		program = 'programs/annovar/', ##annovar has no env program so program has to be reference directly
		datadirectory = 'outputs/annovar_database/',
		version = 'equCab3',
		outName = 'outputs/annovarVCFs/{chrom}/{chrom}_Annovar_SnpEff_Veped',
		outdir = 'outputs/annovarVCFs/{chrom}/',
		unzipped = 'outputs/annovarVCFs/{chrom}/{chrom}_Annovar_SnpEff_Veped.equCab3_multianno.vcf'
	resources:
		time=840,
		mem_mb=100000
	threads: 6
	shell:
		'''
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


##output VCFname needs ot be changed to something generic or a variable
rule combine:
	input: 
		gather_vcfs
	output:
		finalvcf = 'outputs/final/{outputVCF}',
		tbi = 'outputs/final/{outputVCF}.tbi'
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
				IndexFeatureFile -I {output.finalvcf}
		'''

##input VCF change needs to be changed to match with whatever change above is going to happen
##create the annotation table that will serve as basis of performance comparisons
rule prep_table:
	input:
		vcf = expand('outputs/final/{outputVCF}', outputVCF=outputVCF),
		tib = expand('outputs/final/{outputVCF}.tbi',outputVCF=outputVCF)
	output:
		table = 'outputs/final/AnnotationTable.tsv'
	resources:
		time = 420,
		mem_mb = 80000
	shell:
		'''
			gatk VariantsToTable \
				-V {input.vcf} \
				-F CHROM -F POS -F REF -F ALT \
				-F CSQ -F ANN \
				-F Func.ensGene -F Gene.ensGene -F GeneDetail.ensGene \
				-F ExonicFunc.ensGene -F AAChange.ensGene \
				-O {output.table}
		'''


##rule to run Comparison python script
#

rule comparison:
	input:
		table = 'outputs/final/AnnotationTable.tsv',
	output:
		report = 'outputs/final/FinalReport.txt'	
	params:
		program = 'programs/Comparing2.0.py',
		dictionary = 'datafiles/VariantClassDictionary.json',
		VEPPrec = 'datafiles/VEP_precedence.txt',
		EffPrec = 'datafiles/SnpEff_precedence.txt'
	resources:
		time = 360,
		mem_mb = 80000
	shell:
		'''
			python {params.program} {input.table} {params.dictionary} {params.VEPPrec} {params.EffPrec}
		'''

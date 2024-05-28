###Read in table file
##Current attempt: compare impact classifications of variants with one transcript
## compare variant class for variants with one transcript


##Eventually will read in a file that contains the dictionary of variant class matches
##

##creating input arguments for the scripts

import argparse
import json


parser = argparse.ArgumentParser(description='Comparing Annotation simularities between annotation programs')

parser.add_argument('annotation_file',type=str, help='Path to file containing a table of annotations from VEP, SNPEff, and ANNOVAR')
args = parser.parse_args()

table =open(args.annotation_file,'rt') #read in annotation table made by GATK
output = open('outputs/final/FinalReport.txt','wt') #output file
table.readline() #get rid of header
error = open('outputs/final/Errors.txt','wt') #error file
##Files to hold mismatching impacts and classes
mismatch = open('outputs/final/MismatchLines.txt','wt') ##file to hold mismatching info for VEP and SnpEff

V_Amismatch = open('outputs/final/V_AMismatch.txt','wt')
S_Amismatch = open('outputs/final/S_AMismatch.txt','wt')
Somemismatch = open('outputs/final/3AnnotatorComparison.txt','wt')

uneven = open('outputs/final/UniqueTranscripts.txt','wt') ##file to hold records of transcripts with no match
changes = open('outputs/final/ClassChanges.txt','wt')
test = open('outputs/final/testoutput.txt','wt')
partial = open('outputs/final/partialMatches.txt','wt')
linenumber=1 
##impact comparisons
same=0
diff=0
##Class comparisons between VEP and SNPEff
sameClass = 0
diffClass = 0
partClass = 0
##class comparisons between VEP and Annovar
V_Asame =0
V_Adiff =0
##class comparisons between SnpEff and Annovar
S_Asame = 0
S_Adiff = 0
##Class comparisons between all three annotators
AllSame=0
SomeDiff = 0
##Count of annovar partial matches
V_Apart = 0
S_Apart = 0
Eff_part_diff = 0
Eff_part_same = 0
partialmatch2A = 0
totalWhiff = 0

simple =0
singleadditive = 0
simplelines =0
multitranscriptlines =0
multitranscript =0
additivemultitranscript =0

astericks = 0
## dictionary containing the common variant_Class and subclasses that will fall under the common class.
classdict = {
	"splicing_variant":['splicing','splicing_variant','splice_acceptor_variant','splice_donor_variant','splice_region_variant','splice_donor_5th_base_variant','splice_polypyrimidine_tract_variant'],
	"non_coding_transcript_variant":["non_coding_transcript_exon_variant","ncRNA_exonic",'ncRNA','non_coding_transcript_variant','non_coding_exon_variant','ncRNA_intronic'],
	"5_prime_UTR_variant":["5_prime_UTR_variant","UTR5",'5_prime_UTR_truncation','5_prime_UTR_premature_start_codon_gain_variant'],
	"3_prime_UTR_variant":["3_prime_UTR_variant","UTR3",'3_prime_UTR_truncation'],
	"intron_variant":['intronic_variant','intron','intronic','intron_variant','conserved_intron_variant'],
	"upstream_variant":['upstream_gene_variant','upstream'],
	"downstream_variant":["downstream_gene_variant","downstream"],
	"intergenic":['intergenic','intergenic_region','intergenic_variant','conserved_intergenic_variant'],
	"frameshift_variant":["frameshift_variant","frameshift_insertion","frameshift_deletion"],
	'stop_gained':['stop_gained','stopgain'],
	'stop_lost':['stop_lost','stoploss'],
	"inframe_insertion":['inframe_insertion','nonframeshift_insertion','disruptive_inframe_insertion','conservative_inframe_insertion'],
	"inframe_deletion":['inframe_deletion','disruptive_inframe_deletion','nonframeshift_deletion','conservative_inframe_deletion'],
	"nonsynonymous_variant":["missense_variant","nonsynonymous_SNV"],
	"synonymous_variant":["synonymous_variant","synonymous_SNV",'start_retained_variant','stop_retained_variant'],
	"start_lost":['start_lost'],
	"initiator_codon_variant":['initiator_codon_variant'],
	"coding_sequence_variant":["coding_sequence_variant"],
	"transcript_ablation":['transcript_ablation','feature_ablation'],
	"protein_altering_variant":['protein_altering_variant'],
	'exon_loss_variant':['exon_loss_variant'],
	'intragenic_variant':['intragenic_variant'],
	}
mismatch.write('Line\tTranscript\tVEP\tSnpEff\n')
V_Amismatch.write('Line\tVEP\tANNOVAR')
S_Amismatch.write('Line\tSnpEff\tANNOVAR')
Somemismatch.write('Line\tSnpEff\tANNOVAR\tVEP')
differences = {}
V_SClasses={}
V_AClasses={}
S_AClasses={}
partialDict={}
partialVADict={}
partialSADict = {}
bdict = {}
checkcount = 0
## Per variant comparison of annotations
for line in table: ##will loop through each line whic represents a single variant
	VEPDict={} #initating dictionary for lines with multi transcripts
	EffDict={} 
	infoList = line.split('\t') ##splitting columsn of info up.
	linenumber +=1 ##keeping track of what line number is being evaluated

	if '*' in infoList[3]: ##variants with * alleles are not evaluated by annotators so they will be removed
		pass
		error.write('Line number '+str(linenumber)+' removed due to * allele\n')
		astericks +=1
	else: ##evaluating only lines without a *
		VEP = infoList[4] #get VEP annotation
		Eff = infoList[5] #get SNPEff annotation
		annovar = '|'.join(map(str,infoList[6:11])).strip() #grab and reconstruct annovar annotation

		ann_annotation = annovar.split('|')
		ann_area = ann_annotation[0]
		if '3b' in ann_area:
			multiarea= ann_area.split('\\x3b')
			if 'exonic' in multiarea:
				multiarea.append(ann_annotation[3])
			commonlist = []
			for item in multiarea:
				for common,namePos in classdict.items():
					if item in namePos:
						commonName = common
						commonlist.append(commonName)
			ann_class = commonlist
			if ann_area in bdict:
				bdict[ann_area] +=1
			else:
				bdict[ann_area] = 1
		elif ann_area == 'exonic':
			ann_class=ann_annotation[3]
			for common,namePos in classdict.items():
				if ann_class in namePos:
					ann_class = common
		else:
			ann_class=ann_area
			for common,namePos in classdict.items():
				if ann_class in namePos:
					ann_class = common
	
		######At this point we start evaluating differently if the variant overlaps multiple transcripts	

		if ',' in VEP and ',' in Eff:  ##finding lines that have multiple transcripts
			TranSplitVEP = VEP.split(',') ##separating individual transcript predictions for VEP and SnpEff
			TranSplitEff = Eff.split(',')

			multitranscriptlines += 1
			##Checking if there are different numbers of transcripts being used. SNPeff will typically have more	
			if len(TranSplitVEP) != len(TranSplitEff): ##making note of lines with odd numbers of transcripts
				error.write('Line number ' +str(linenumber) + ' has a mismatched number of transcripts\n'
				+ 'VEP has '+str(len(TranSplitVEP))+ ' transcripts and SnpEff has '+str(len(TranSplitEff))+'\n\n')

			#if '3b' in ann_class:
			#	multinoadditive3b +=1
	
			##creating dictionaries for the multiple transcripts of VEP and SNPEff. Wills save transcript:class:impact
			for entryNum in range(len(TranSplitVEP)): #create dictionary for VEP transcripts
				CurrentVEP = TranSplitVEP[entryNum] ##grabs current transcript/annotation from list of transcript/annotations using index of list
				VEP_annotation = CurrentVEP.split('|') ##split annotation report into individual info fields
				VEPDict[VEP_annotation[6]]={'class':VEP_annotation[1],'impact':VEP_annotation[2]} # create dict entry to pair transcript with variant class and impact prediction
			for entryNum in range(len(TranSplitEff)): ##mirror of above but with snpeff
				CurrentEff = TranSplitEff[entryNum]
				Eff_annotation = CurrentEff.split('|')
				EffDict[Eff_annotation[6].split('.')[0]]={'class':Eff_annotation[1],'impact':Eff_annotation[2]} #transcript number contains version number that is removed for pairing


			##Determining if ther are transcripts not shared across VEP and SNPEff
			VEPunique = set(VEPDict.keys()) - set(EffDict.keys()) #Determining transcripts that one tool uses that the other tool does not
			Effunique = set(EffDict.keys()) - set(VEPDict.keys())

			##create reports of any nonshared transcripts
			if len(VEPunique) > 0: #creating records of transcripts that have no match
				for VU in VEPunique:
					print("VEP has a unique transcript:",VU,"with class", VEPDict[VU]['class'], "and impact",VEPDict[VU]['impact'],".\n",file=uneven)
			if len(Effunique) > 0: #creating records of transcripts that have no match
				for EU in Effunique:
					print("SnpEff has a unique transcript:",EU,"with class", EffDict[EU]['class'], "and impact",EffDict[EU]['impact'],".\n",file=uneven)
			

			#### Making comparisons of annotations from transcripts that are found in both VEP and SnpEff
			SharedTranscripts = set(VEPDict.keys()) & set(EffDict.keys()) ##creating list of only transcripts present in both lists
			for transcript in SharedTranscripts:


				###Some annotations will have additive annotations. '&'
				###Annotations with & have to be split specially.
				##This top section will only be done on annotations with no &
				if '&' not in VEPDict[transcript]['class'] and '&' not in EffDict[transcript]['class']: ##passes on annotations that have & in either VEP or SNPEff

					
					multitranscript +=1
					#must normalize variant classification names between annotator programs
					## From both dictionaries grab the impact and class predictions made with associated transcript ID
					for commonName, NamePos in classdict.items(): #iterating through key:value pairs in dictionary
						MatchFound = False
						if VEPDict[transcript]['class'] in NamePos: ##if the classification from VEP is in the list of possibilities for the current common name
							VEPDict[transcript]['common'] = commonName ##record common name in dictionary for this transcript
							MatchFound = True ##set match to true
							break ##stop cycle

					#this will check for any classes that don't have common names. This should not occur anymore
					if MatchFound == False:
						VEPDict[transcript]['common']=VEPDict[transcript]['class']
						print('No match found for VEP variant class',VEPDict[transcript]['class'],'Maintaing original classification',file=error)


					###Same as above but for SNPEff
					for commonName, NamePos in classdict.items():
						MatchFound = False
						if EffDict[transcript]['class'] in NamePos:
							EffDict[transcript]['common'] = commonName
							MatchFound = True
							break
					if MatchFound == False: 
						EffDict[transcript]['common']=EffDict[transcript]['class']
						print('No match found for SnpEff variant class', EffDict[transcript]['class'],'Maintaining original classification',file=error)
					
					##Following normalization record the impact and class corresponding to each tool of the current transcript
					VEP_impact = VEPDict[transcript]['impact']
					VEP_class = VEPDict[transcript]['common']
					Eff_impact = EffDict[transcript]['impact'] 
					Eff_class = EffDict[transcript]['common']
					#print('This is the impact pair '+VEP_impact+' = '+Eff_impact+'\n')
					

					#Comparing Impact predictions and increase count of match or mismatched
					if VEP_impact == Eff_impact:
						#print('Impacts Match!')
						same += 1 ##increases count of the amount of identical impact predictions
					else:
						#print('They do not match!')
			
						##print out mismatching records for checking later
						#print(str(linenumber),transcript,VEP_impact,Eff_impact,sep='\t',file=mismatch) ##printing the mismatches
						diff += 1 ## increases count of the number of disagreeing impact predictions
						difference = VEP_impact+'_'+Eff_impact ##keeping track of which kind of mismatch is occuring LOW_MODERATE etc
						if difference in differences: ## if we've already seen one increase count of it
							differences[difference] += 1
						else:
							differences[difference] = 1 ##if not make new categor
					
					#Comparing variant class assignments. Increase count of match or mistmatch

	
					###Code to keep track of what variant_classes and how many are found
					#if '&' in VEP_class:
					#	ClassList = VEP_class.split('&')
					#	for ann in ClassList:
					#		if ann in Classes:
					#			Classes[ann] += 1
					#		else:
					#			Classes[ann] =1
					#else:
					#	if VEP_class in Classes:
					#		Classes[VEP_class] +=1
					#	else:
					#		Classes[VEP_class]=1


					#if '&' in Eff_class:
					#	ClassList = Eff_class.split('&')
					#	for ann in ClassList:
					#		if ann in Classes:
					#			Classes[ann] += 1
					#		else:	
					#			Classes[ann] =1
					#else:
					#	if Eff_class in Classes:
					#		Classes[Eff_class] +=1
					#	else:
					#		Classes[Eff_class]=1
					

					##comparing variant classes. Again these do not have '&'
					if VEP_class == Eff_class:
						#print('The classes match!')
						sameClass += 1
					else:
						#print('The classes do not match!')
						diffClass += 1
						hold = VEP_class+'_'+Eff_class
						if hold in V_SClasses:
							V_SClasses[hold] +=1
						else:
							V_SClasses[hold] = 1
						#print(str(linenumber),transcript,VEP_class,Eff_class,sep='\t',file=mismatch)

					if isinstance(ann_class,list):
						#print('length of annovar is more than 1',ann_class,file=test)
						if VEP_class in ann_class:
							V_Asame+=1
						else:
							V_Adiff +=1
							#print("VEP",VEP_class,ann_class, file =test)
						if Eff_class in ann_class:
							S_Asame +=1
						else:
							S_Adiff+=1
							#print("VEP",Eff_class,ann_class,file=test)
						checkcount +=1
					else:
					#	print('length of annovar is equal to 1', ann_class, file = test)
						if VEP_class == ann_class:
							V_Asame +=1
						else:
							V_Adiff +=1
							hold = VEP_class+'_'+ann_class
							if hold in V_AClasses:
								V_AClasses[hold] +=1
							else:
								V_AClasses[hold] = 1

						if Eff_class == ann_class:
							S_Asame +=1
						else:
							S_Adiff +=1
							hold = Eff_class+'_'+ann_class
							if hold in S_AClasses:
								S_AClasses[hold] +=1
							else:
								S_AClasses[hold] = 1

					
					if Eff_class in ann_class and VEP_class in ann_class:
						AllSame += 1
					else:
						SomeDiff +=1
					
				else: ###for annotations that are additive and have a &
					##Quickly compare effects in the same method as before
					additivemultitranscript += 1
					VEP_impact = VEPDict[transcript]['impact']
					Eff_impact = EffDict[transcript]['impact']
					if '3b' in ann_class:
						multiadditive3b += 1
					if VEP_impact == Eff_impact:
						same += 1
					else:
						diff += 1
						#print(str(linenumber),transcript,VEP_impact,Eff_impact,sep='\t',file=mismatch)
						difference = VEP_impact+'_'+Eff_impact
						if difference in differences:
							differences[difference] += 1
						else:
							differences[difference] = 1

					##Gather class information
					VEP_class = VEPDict[transcript]['class']
					Eff_class = EffDict[transcript]['class']
						
					##split the annotations on the & and create lists with all possibilities
					Vposs = VEP_class.split('&')
					Eposs = Eff_class.split('&')

					##Include loop to change all items in lists to common names
					commonVposs = []
					commonEposs = []
					for item in Vposs: #for all the annotations in VEP
						for common,namePos in classdict.items():
							if item in namePos: ##if the annotation is within a common category
								commonVposs.append(common) ##add common name to new list
					##add list to dictionary under 'common'
					VEPDict[transcript]['common'] = commonVposs

					for item in Eposs: ##same as above
						for common,namePos in classdict.items():
							if item in namePos:
								commonEposs.append(common)
					
					##add common name list to dictionary under 'common'
					EffDict[transcript]['common'] = commonEposs
					##find overlap of both lists of annotations. These are matches
					overlap = list(set(commonVposs) & set(commonEposs))
		
					##find values only found in one list of the other. These are mismatches
					Veponly = [x for x in commonVposs if x not in commonEposs] ##for every item in VEP if this item is not found in Eff add to Veponly
					Effonly = [x for x in commonEposs if x not in commonVposs]
					unique = len(Veponly)+len(Effonly) ##determining total number of classes unique to either list
					match = 100*(len(overlap)/(unique+len(overlap))) ##determine the percent of annotations that match
					if match == 100: ##if the match is 100 it will be considered a full match and added to same class count
						sameClass +=1
					elif match == 0:
						diffClass += 1
						#hold = str(Veponly)+'_'+str(Effonly)
						#if hold in V_SClasses:
						#	V_SClasses[hold] += 1
						#else:
						#	V_SClasses[hold] =1
					else:
						partClass += 1 ##if it does not match it will be added to the diff class count AND...

						##print the partial matches for records
						print(match,file=partial)
						print(overlap,Veponly,Effonly,sep='\t',file=partial)
						##keep track of what different levels of partial matches are being made
						if match in partialDict:
							partialDict[match] +=1
						else:
							partialDict[match] =1
						#hold = str(Veponly)+'_'+str(Effonly)
						#if hold in V_SClasses:
						#	V_SClasses[hold] += 1
						#else:
						#	V_SClasses[hold] =1
					####This will be changed in the future for multiple transcripts of ANNOVAR
					if isinstance(ann_class,str):
						if ann_class in VEPDict[transcript]['common']:
							length = len(VEPDict[transcript]['common'])
							match = 100*(1/length)
							if match ==100:
								V_Asame += 1
							else:	
								V_Apart +=1
								if match in partialVADict:
									partialVADict[match] +=1
								else:
									partialVADict[match] =1
						else:
							
							V_Adiff +=1
						

						if ann_class in EffDict[transcript]['common']:
							length = len(EffDict[transcript]['common'])
							match = 100*(1/length)
							if match ==100:
								S_Asame +=1
							else:
								S_Apart +=1
								if match in partialSADict:
									partialSADict[match] += 1
								else:
									partialSADict[match] =1
						else:
							
							S_Adiff +=1
					else:
						VAoverlap = list(set(commonVposs) & set(ann_class))
						SAoverlap = list(set(commonEposs) & set(ann_class))
						
						VEPonly = [x for x in commonVposs if x not in ann_class]
						annonly = [x for x in ann_class if x not in commonVposs]
						unique = len(VEPonly)+len(annonly)
						match = 100*(len(VAoverlap)/(unique+len(VAoverlap)))

						if match == 100:
							V_Asame +=1
						elif match ==0:
							V_Adiff +=1
						else:
							V_Apart +=1
						Effonly = [x for x in commonEposs if x not in ann_class]
						annonly = [x for x in ann_class if x not in commonEposs]
						unique = len(Effonly)+len(annonly)
						match = 100*(len(SAoverlap)/(unique+len(SAoverlap)))
						if match ==100:
							S_Asame +=1
						elif match ==0:
							S_Adiff +=1
						else:
							S_Apart +=1


			#if '3b' in ann_class:
				#print(ann_class,file=test)
				#print(line,file=test)
				#if ann_class in bdict:
				#	bdict[ann_class] +=1
				#else:
				#	bdict[ann_class] =1


		##COMPARISONS BEING MADE FOR VARIANTS WITH ONLY ONE TRANSCRIPT
		else:
			simplelines +=1
			#output.write(VEP+'\t'+Eff+'\t'+annovar+'\n')
			##Since there's only one transcript can jump straight into is
			VEP_annotation = VEP.split('|') #split VEP annotation into indv. info fields
			Eff_annotation = Eff.split('|') #split snpEff annotations into inv. info fields
			#annovar_annotation = annovar.split('|')
			VEP_transcript= VEP_annotation[6]
			EFF_transcript= Eff_annotation[6].split('.')[0]

			##error check. This will be true often for intergenic regions
			if VEP_transcript != EFF_transcript:
				print('Error: Transcripts do not match at line', linenumber,VEP_transcript,EFF_transcript,file=error)

			VEP_class = VEP_annotation[1]
			VEP_impact = VEP_annotation[2]
			Eff_class = Eff_annotation[1]
			Eff_impact = Eff_annotation[2]
			#ann_annotation = annovar.split('|')
			#ann_area = ann_annotation[0]
			

			##annovar has a field for area and then an extra area for exonic function if the area is exonic. This will make sure the correct class is grabbed
			#if ann_area == 'exonic':
			#	ann_class = ann_annotation[3]
			#else:
			#	ann_class = ann_area


			##Changing the annovar class to the common name
			#for common,namePos in classdict.items():
			#	if ann_class in namePos:
			#		ann_class = common	



			##Comparing differently if there is & or no &
			if '&' in VEP_class or '&' in Eff_class: ##If there are mutliple annotations will need to split and perform partial matches

				singleadditive +=1
				###split
				Vposs = VEP_class.split('&')
				Eposs = Eff_class.split('&')
				commonVposs = []
				commonEposs = []
				if '3b' in ann_class:
					singleadditive3b +=1
				### iterate and change names
				for item in Vposs:
					for common,namePos in classdict.items():
						if item in namePos:
							commonVposs.append(common)
				for item in Eposs:
					for common,namePos in classdict.items():
						if item in namePos:
							commonEposs.append(common)


				###perform partial match checks
				overlap = list(set(commonVposs) & set(commonEposs))
				Veponly = [x for x in commonVposs if x not in commonEposs]
				Effonly = [x for x in commonEposs if x not in commonVposs]
				unique = len(Veponly)+len(Effonly)
				match = 100*(len(overlap)/(unique+len(overlap)))
				if match == 100:
					sameClass += 1
				elif match == 0:
					diffClass += 1
				else:
					partClass +=1
					print(match,file=partial)
					print(overlap,Veponly,Effonly,sep='\t',file=partial)
					if match in partialDict:
						partialDict[match] +=1
					else:
						partialDict[match] = 1
				
				if isinstance(ann_class,str):
			
					if ann_class in commonVposs:
						match = 100*(1/len(commonVposs))
						if match ==100:
							V_Asame +=1
						else:
							V_Apart +=1
							if match in partialVADict:
								partialVADict[match] +=1
							else:
								partialVADict[match]=1
					else:
						V_Adiff +=1

					if ann_class in commonEposs:
						match = 100*(1/len(commonEposs))
						if match ==100:
							S_Asame +=1
						else:
							S_Apart +=1
							if match in partialSADict:
								partialSADict[match]+=1
							else:
								partialSADict[match]=1
					else:
						S_Adiff+=1
				else:
					VAoverlap = list(set(commonVposs) & set(ann_class))
					SAoverlap = list(set(commonEposs) & set(ann_class))

					VEPonly = [x for x in commonVposs if x not in ann_class]
					annonly = [x for x in ann_class if x not in commonVposs]
					unique = len(VEPonly)+len(annonly)
					match = 100*(len(VAoverlap)/(unique+len(VAoverlap)))

					if match ==100:
						V_Asame +=1
					elif match ==0:
						V_Adiff +=1
					else:
						V_Apart +=1
					Effonly = [x for x in commonEposs if x not in ann_class]
					annonly = [x for x in ann_class if x not in commonEposs]
					unique = len(Effonly)+len(annonly)
					match = 100*(len(SAoverlap)/(unique+len(SAoverlap)))
					if match ==100:
						S_Asame +=1
					elif match ==0:
						S_Adiff += 1
					else:
						S_Apart +=1
					
				## I also keep track of how many times this 'partial total' match is needed
				##might be a good idea to eventually code which other classes it is not matching with
				#if 'x3b' in ann_class:
				#	print('x3b in yes & one transcript',file=test)
				#	print(ann_class,line,file=test)
				#if ann_class in commonVposs:
				#	V_Asame += 1
				#	#print("Same",ann_class, commonVposs,file=test)
				#	VEP_part_same += 1
				#else:
				#	V_Adiff += 1
				#	#print("Different",ann_class, commonVposs,file=test)
				#	VEP_part_diff += 1
				#if ann_class in commonEposs:
				#	S_Asame += 1
				#	Eff_part_same += 1
				#else:
				#	S_Adiff += 1
				#	Eff_part_diff += 1
				#if ann_class in commonVposs and ann_class in commonEposs:
				#	AllSame += 1
				#	partialmatch2A +=1
				#else:
				#	SomeDiff += 1
				#	totalWhiff += 1

					


			else: ##if there aren't any & can proceed as normal
				##Changing the names of all the variant classes for all programs
				if '3b' in ann_class:
					singlenoadditive3b +=1
				simple += 1
				for commonName, NamePos in classdict.items():
					if VEP_annotation[1] in NamePos:
						VEP_class = commonName
						#print('VEP classification',VEP_annotation[1],'has been changed to',VEP_class,file=changes)
					if Eff_annotation[1] in NamePos:
						Eff_class = commonName
						#print('SnpEff classification',Eff_annotation[1],'has been changed to',Eff_class,file=changes)
					if ann_class in NamePos:
						#print('Annovar classification',ann_class,'has been changed to',commonName,file=changes)
						ann_class =commonName
				#output.write(VEP_class+'\t'+VEP_impact+'\t'+Eff_class+'\t'+Eff_impact+'\t'+annovar_class+'\n')

				
				###Code to keep track of what variant classes and how many are being found. This will have to be moved as it is currently in a section where there are no &s passing
#				if '&' in VEP_class:
#					ClassList = VEP_class.split('&')
#					for ann in ClassList:
#						if ann in Classes:
#							Classes[ann] += 1
#						else:
#							Classes[ann] =1
#				else:
#					if VEP_class in Classes:
#						Classes[VEP_class] +=1
#					else:
#						Classes[VEP_class]=1
#				if '&' in Eff_class:
#					ClassList = Eff_class.split('&')
#					for ann in ClassList:
#						if ann in Classes:
#							Classes[ann] += 1
#						else:
#							Classes[ann] =1
#				else:
#					if Eff_class in Classes:
#						Classes[Eff_class] +=1
#					else:
#						Classes[Eff_class]=1


				#if 'x3b' in ann_class:
				#	print('no & one transcript',file=test)
				#	print(line,file=test)
				#	if ann_class in bdict:
				#		bdict['ann_class'] +=1
				#	else:
				#		bdict['ann_class']=1
					
				##Comparing VEP and Eff class
				if VEP_class == Eff_class:
					sameClass +=1
				else:
					diffClass +=1
					hold = VEP_class+'_'+Eff_class
					if hold in V_SClasses:
						V_SClasses[hold] +=1
					else:
						V_SClasses[hold]=1
					#print(str(linenumber),'NA',VEP_class,Eff_class,sep='\t',file=mismatch)
			
				##Comparing VEP and Ann
				if VEP_class == ann_class:
					V_Asame +=1
				else:
					V_Adiff +=1
					print('VEP',VEP_class,ann_class,file=test)
					hold = VEP_class+'_'+ann_class
					if hold in V_AClasses:
						V_AClasses[hold]+=1
					else:
						V_AClasses[hold]=1
					print(linenumber,VEP_class,ann_class,sep='\t',file=V_Amismatch)

				##Comparing Eff and Ann
				if Eff_class == ann_class:
					S_Asame +=1
				else:
					S_Adiff +=1
					print('Eff',Eff_class,ann_class,file=test)
					hold = Eff_class+'_'+ann_class
					if hold in S_AClasses:
						S_AClasses[hold] +=1
					else:
						S_AClasses[hold]=1
					print(linenumber,Eff_class,ann_class,sep='\t',file=S_Amismatch)


				##Comparing all three
				if Eff_class == ann_class == VEP_class:
					AllSame +=1
				else:
					SomeDiff +=1
					print(linenumber,Eff_class,ann_class,VEP_class,sep='\t',file=Somemismatch)
			##impact
				
			if VEP_impact == Eff_impact:
				same +=1
			else: 
				diff +=1
				#print(str(linenumber),VEP_transcript,VEP_impact,Eff_impact,sep='\t',file=mismatch)
				difference = VEP_impact+'_'+Eff_impact
				if difference in differences:
					differences[difference] +=1 
				else:
					differences[difference]=1





####Creating report file
ImpactMatchPer = 100*(same/(same+diff))
ImpactMisMatchPer = 100*(diff/(same+diff))
ClassMatchPer = 100*(sameClass/(sameClass+diffClass+partClass))
ClassMisMatchPer = 100*(diffClass/(sameClass+diffClass+partClass))
ClassPartMatch = 100*(partClass/(sameClass+diffClass+partClass))

ClassMatchV_A = 100*(V_Asame/(V_Asame+V_Adiff+V_Apart))
ClassMatchS_A = 100*(S_Asame/(S_Asame+S_Adiff+S_Apart))
ClassMatchAll = 100*(AllSame/(AllSame+SomeDiff))
ClassPartMatchS_A = 100*(S_Apart/(S_Asame+S_Adiff+S_Apart))
ClassPartMatchV_A = 100*(V_Apart/(V_Asame+V_Adiff+V_Apart))

ClassMismatchV_A = 100*(V_Adiff/(V_Asame+V_Adiff))
ClassMismatchS_A = 100*(S_Adiff/(S_Asame+S_Adiff))
ClassMismatchSome = 100*(SomeDiff/(AllSame+SomeDiff))


output.write('There were ' +str(same)+ ' matching impact predictions ('+ str(ImpactMatchPer)+ '%).\nThere were ' + str(diff) + ' different impact predictions ('+str(ImpactMisMatchPer)+'%).\n')
output.write('VEP_SnpEff\nThere were ' +str(sameClass)+ ' matching variant class predictions ('+str(ClassMatchPer)+ '%).\nThere were '+str(partClass)+ ' partial variant class predictions ('+str(ClassPartMatch)+'\n'+'There were ' + str(diffClass) + ' different variant class predictions ('+ str(ClassMisMatchPer)+ '%).\n')


output.write('VEP_ANNOVAR\nThere were '+str(V_Asame)+ ' matching variant class predictions ('+str(ClassMatchV_A)+'%).\n'+'There were '+str(V_Apart)+' partial variant class predictions ('+str(ClassPartMatchV_A)+'%).\n'+'There were ' + str(V_Adiff) + 'difference variant class predictions ('+str(ClassMismatchV_A)+'%).\n')


output.write('SnpEff_ANNOVAR\nThere were '+str(S_Asame)+ ' matching variant class predictions ('+str(ClassMatchS_A)+'%).\n'+'There were'+str(S_Apart)+' partial variant class predictions ('+str(ClassPartMatchS_A)+'%).\n'+'There were ' + str(S_Adiff) + 'different variant class predictions ('+str(ClassMismatchS_A)+'%).\n')


output.write('VEP_ANNOVAR_SnpEff\nThere were '+str(AllSame)+ ' matching variant class predictions ('+str(ClassMatchAll)+'%).\nThere were ' + str(SomeDiff) + 'different variant class predictions ('+str(ClassMismatchSome)+'%).\n')

print('VEP Impact\tSnpEff Impact\tNumber of Occurences',file=output)
for key in differences.keys():
	first,second = key.split('_')
	occurence = differences[key]
	print(first,second,occurence,sep='\t',file=output)
#for key in Classes.keys():
#	print(key,Classes[key],file=output)
print("Partial States VEP_SNPEff",file=output)
for key in partialDict.keys():
	print(key,partialDict[key],file=output)
print("Partial Stats VEP_Annovar",file=output)
for key in partialVADict.keys(): 
        print(key,partialVADict[key],file=output)
print("Partial Stats SNPEff_Annovar",file=output)
for key in partialSADict.keys(): 
	print(key,partialSADict[key],file=output)
for key in bdict.keys():
	print(key,bdict[key],file=output)


print('VEP','SNPEff','Count',sep='\t',file=mismatch)
for key in V_SClasses.keys():
	print(key,V_SClasses[key],sep='\t',file=mismatch)
print('\n\n\n\n',file=mismatch)
print('VEP','ANNOVAR','Count',sep='\t',file=mismatch)
for key in V_AClasses.keys():
        print(key,V_AClasses[key],sep='\t',file=mismatch)
print('\n\n\n\n',file=mismatch)
print('SNPEff','ANNOVAR','Count',sep='\t',file=mismatch)
for key in S_AClasses.keys():
        print(key,S_AClasses[key],sep='\t',file=mismatch)


#print('Total Difference from VEP:',VEP_part_diff,"Partial Match with VEP:",VEP_part_same,"Total Difference from Eff:",Eff_part_diff,"Partial Match with Eff:",Eff_part_same,"Matched Nothing:",totalWhiff,sep='\n',file=test)
print('Total lines with multiple transcripts',multitranscriptlines,'Comparisons from multiple transcript lines with annotations that are simple',multitranscript,'Multitranscript additive comparisons',additivemultitranscript,'Total lines with single transcripts',simplelines,'Single transcript non additive comparisons',simple,'Single transcript with additive annotations',singleadditive,sep='\n',file=output)
print(checkcount,file=output)

table.close()
output.close()
error.close()
mismatch.close()
uneven.close()
changes.close()
test.close()
V_Amismatch.close()
S_Amismatch.close()
Somemismatch.close()
partial.close()

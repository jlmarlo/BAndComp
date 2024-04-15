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

#with open(ClassDictionary.json) as f:
#	classdict = json.load(f)
table =open(args.annotation_file,'rt') #read in annotation table made by GATK
output = open('FinalReport.txt','wt') #output file
table.readline() #get rid of header
error = open('Errors.txt','wt') #error file
##Files to hold mismatching impacts and classes
mismatch = open('MismatchLines.txt','wt') ##file to hold mismatching info for VEP and SnpEff

V_Amismatch = open('V_AMismatch.txt','wt')
S_Amismatch = open('S_AMismatch.txt','wt')
Somemismatch = open('3AnnotatorComparison.txt','wt')

uneven = open('UniqueTranscripts.txt','wt') ##file to hold records of transcripts with no match
changes = open('ClassChanges.txt','wt')
test = open('testoutput.txt','wt')
linenumber=1 
##impact comparisons
same=0
diff=0
##Class comparisons between VEP and SNPEff
sameClass = 0
diffClass = 0
##class comparisons between VEP and Annovar
V_Asame =0
V_Adiff =0
##class comparisons between SnpEff and Annovar
S_Asame = 0
S_Adiff = 0
##Class comparisons between all three annotators
AllSame=0
SomeDiff = 0


classdict = {
	"intergenic":['intergenic','intergenic_region','intergenic_variant'],
	"inframe_deletion":['inframe_deletion','disruptive_inframe_deletion'],
	"intronic_variant":['intronic_variant','intron','intronic','intron_variant']}
mismatch.write('Line\tTranscript\tVEP\tSnpEff\n')
V_Amismatch.write('Line\tVEP\tANNOVAR')
S_Amismatch.write('Line\tSnpEff\tANNOVAR')
Somemismatch.write('Line\tSnpEff\tANNOVAR\tVEP')
differences = {}
## Per variant comparison of annotations
for line in table:
	VEPDict={} #initating dictionary for lines with multi transcripts
	EffDict={} 
	infoList = line.split('\t') ##splitting columsn of info up
	linenumber +=1 ##keeping track of what line number is being evaluated

	if '*' in infoList[3]: ##variants with * alleles are not evaluated by annotators so they will be removed
		pass
		error.write('Line number '+str(linenumber)+' removed due to * allele\n')

	else:
		VEP = infoList[4] #get VEP annotation
		Eff = infoList[5] #get SNPEff annotation
		annovar = '|'.join(map(str,infoList[6:11])).strip() #grab and reconstruct annovar annotation
		print(annovar,file=test)


##Code for comparisons between VEP and SNPEff



		if ',' in VEP and ',' in Eff:  ##finding lines that have multiple transcripts
			TranSplitVEP = VEP.split(',') ##separating individual transcript predictions
			TranSplitEff = Eff.split(',')
			if len(TranSplitVEP) != len(TranSplitEff): ##making note of lines with odd numbers of transcripts
				error.write('Line number ' +str(linenumber) + ' has a mismatched number of transcripts\n'
				+ 'VEP has '+str(len(TranSplitVEP))+ ' transcripts and SnpEff has '+str(len(TranSplitEff))+'\n\n')
			for entryNum in range(len(TranSplitVEP)): #create dictionary for VEP transcripts
				CurrentVEP = TranSplitVEP[entryNum]
				VEP_annotation = CurrentVEP.split('|') ##split annotation report into individual info fields
				VEPDict[VEP_annotation[6]]={'class':VEP_annotation[1],'impact':VEP_annotation[2]} # create dict entry to pair transcript with variant class and impact prediction
			for entryNum in range(len(TranSplitEff)):
				CurrentEff = TranSplitEff[entryNum]
				Eff_annotation = CurrentEff.split('|')
				EffDict[Eff_annotation[6].split('.')[0]]={'class':Eff_annotation[1],'impact':Eff_annotation[2]} #transcript number contains version number that is removed for pairing
			VEPunique = set(VEPDict.keys()) - set(EffDict.keys()) #Determining transcripts that one tool uses that the other tool does not
			Effunique = set(EffDict.keys()) - set(VEPDict.keys())
			if len(VEPunique) > 0: #creating records of transcripts that have no match
				for VU in VEPunique:
					print("VEP has a unique transcript:",VU,"with class", VEPDict[VU]['class'], "and impact",VEPDict[VU]['impact'],".\n",file=uneven)
			if len(Effunique) > 0: #creating records of transcripts that have no match
				for EU in Effunique:
					print("SnpEff has a unique transcript:",EU,"with class", EffDict[EU]['class'], "and impact",EffDict[EU]['impact'],".\n",file=uneven)
			


			SharedTranscripts = set(VEPDict.keys()) & set(EffDict.keys()) ##creating list of only transcripts present in both lists
			for transcript in SharedTranscripts:
				## From both dictionaries grab the impact and class predictions made with associated transcript ID
				for commonName, NamePos in classdict.items():
					MatchFound = False
					if VEPDict[transcript]['class'] in NamePos:
						VEPDict[transcript]['common'] = commonName
						print('VEP classification',VEPDict[transcript]['class'], 'has been changed to',VEPDict[transcript]['common'],file = changes)
						MatchFound = True
						break
				if MatchFound == False:
					VEPDict[transcript]['common']=VEPDict[transcript]['class']
					print('No match found for VEP variant class',VEPDict[transcript]['class'],'Maintaing original classification',file=changes)
				for commonName, NamePos in classdict.items():
					MatchFound = False
					if EffDict[transcript]['class'] in NamePos:
						EffDict[transcript]['common'] = commonName
						print('SnpEff classification',EffDict[transcript]['class'],'has been changed to',EffDict[transcript]['common'],file=changes)
						MatchFound = True
						break
				if MatchFound == False: 
					EffDict[transcript]['common']=EffDict[transcript]['class']
					print('No match found for SnpEff variant class', EffDict[transcript]['class'],'Maintaining original classification',file=changes)

				VEP_impact = VEPDict[transcript]['impact']
				VEP_class = VEPDict[transcript]['common']
				Eff_impact = EffDict[transcript]['impact'] 
				Eff_class = EffDict[transcript]['common']
				#print('This is the impact pair '+VEP_impact+' = '+Eff_impact+'\n')
				#Comparing Impact predictions and increase count of match or mismatched
				if VEP_impact == Eff_impact:
					#print('Impacts Match!')
					same += 1
				else:
					#print('They do not match!')
					##print out mismatching records for checking later
					print(str(linenumber),transcript,VEP_impact,Eff_impact,sep='\t',file=mismatch)
					diff += 1
					difference = VEP_impact+'_'+Eff_impact
					if difference in differences:
						differences[difference] += 1
					else:
						differences[difference] = 1
				#Comparing variant class assignments. Increase count of match or mistmatch	
				if VEP_class == Eff_class:
					#print('The classes match!')
					sameClass += 1
				else:
					#print('The classes do not match!')
					diffClass += 1
					print(str(linenumber),transcript,VEP_class,Eff_class,sep='\t',file=mismatch)
				#print('This is the class pair ' +VEP_class+' = '+Eff_class+'\n')				
			###This will be where I incorporate ANNOVAR comparison when VEP and SNPEff have multiple transcripts
			##I will have to find a way to pull the most deleterious variant class from the VEP and SNPEff dictionaries
			##And then compare that to the one that annovar has





		else:
			#output.write(VEP+'\t'+Eff+'\t'+annovar+'\n')
			VEP_annotation = VEP.split('|') #split VEP annotation into indv. info fields
			Eff_annotation = Eff.split('|') #split snpEff annotations into inv. info fields
			annovar_annotation = annovar.split('|')
			VEP_transcript= VEP_annotation[6]
			EFF_transcript= Eff_annotation[6].split('.')[0]
			if VEP_transcript != EFF_transcript:
				print('Error: Transcripts do not match at line', linenumber,VEP_transcript,EFF_transcript,file=error)

			VEP_class = VEP_annotation[1]
			VEP_impact = VEP_annotation[2]
			Eff_class = Eff_annotation[1]
			Eff_impact = Eff_annotation[2]
			ann_annotation = annovar.split('|')
			ann_area = ann_annotation[0]
			if ann_area == 'exonic':
				ann_class = ann_annotation[3]
			else:
				ann_class = ann_area
			
			for commonName, NamePos in classdict.items():
				if VEP_annotation[1] in NamePos:
					VEP_class = commonName
					print('VEP classification',VEP_annotation[1],'has been changed to',VEP_class,file=changes)
				if Eff_annotation[1] in NamePos:
					Eff_class = commonName
					print('SnpEff classification',Eff_annotation[1],'has been changed to',Eff_class,file=changes)
				if ann_class in NamePos:
					print('Annovar classification',ann_class,'has been changed to',commonName,file=changes)
					ann_class =commonName
			#output.write(VEP_class+'\t'+VEP_impact+'\t'+Eff_class+'\t'+Eff_impact+'\t'+annovar_class+'\n')
			if VEP_impact == Eff_impact:
				same +=1
			else: 
				diff +=1
				print(str(linenumber),VEP_transcript,VEP_impact,Eff_impact,sep='\t',file=mismatch)
				difference = VEP_impact+'_'+Eff_impact
				if difference in differences:
					differences[difference] +=1 
				else:
					differences[difference]=1
			if VEP_class == Eff_class:
				sameClass +=1
			else:
				diffClass +=1
				print(str(linenumber),'NA',VEP_class,Eff_class,sep='\t',file=mismatch)
			if VEP_class == ann_class:
				V_Asame +=1
			else:
				V_Adiff +=1
				print(linenumber,VEP_class,ann_class,sep='\t',file=V_Amismatch)
			if Eff_class == ann_class:
				S_Asame +=1
			else:
				S_Adiff +=1
				print(linenumber,Eff_class,ann_class,sep='\t',file=S_Amismatch)
			if Eff_class == ann_class == VEP_class:
				AllSame +=1
			else:
				SomeDiff +=1
				print(linenumber,Eff_class,ann_class,VEP_class,sep='\t',file=Somemismatch)


####Creating report file
ImpactMatchPer = 100*(same/(same+diff))
ImpactMisMatchPer = 100*(diff/(same+diff))
ClassMatchPer = 100*(sameClass/(sameClass+diffClass))
ClassMisMatchPer = 100*(diffClass/(sameClass+diffClass))

ClassMatchV_A = 100*(V_Asame/(V_Asame+V_Adiff))
ClassMatchS_A = 100*(S_Asame/(S_Asame+S_Adiff))
ClassMatchAll = 100*(AllSame/(AllSame+SomeDiff))

ClassMismatchV_A = 100*(V_Adiff/(V_Asame+V_Adiff))
ClassMismatchS_A = 100*(S_Adiff/(S_Asame+S_Adiff))
ClassMismatchSome = 100*(SomeDiff/(AllSame+SomeDiff))


output.write('There were ' +str(same)+ ' matching impact predictions ('+ str(ImpactMatchPer)+ '%).\nThere were ' + str(diff) + ' different impact predictions ('+str(ImpactMisMatchPer)+'%).\n')
output.write('VEP_SnpEff\nThere were ' +str(sameClass)+ ' matching variant class predictions ('+str(ClassMatchPer)+ '%).\nThere were ' + str(diffClass) + ' different variant class predictions ('+ str(ClassMisMatchPer)+ '%).\n')
output.write('VEP_ANNOVAR\nThere were '+str(V_Asame)+ ' matching variant class predictions ('+str(ClassMatchV_A)+'%).\nThere were ' + str(V_Adiff) + 'difference variant class predictions ('+str(ClassMismatchV_A)+'%).\n')
output.write('SnpEff_ANNOVAR\nThere were '+str(S_Asame)+ ' matching variant class predictions ('+str(ClassMatchS_A)+'%).\nThere were ' + str(S_Adiff) + 'different variant class predictions ('+str(ClassMismatchS_A)+'%).\n')
output.write('VEP_ANNOVAR_SnpEff\nThere were '+str(AllSame)+ ' matching variant class predictions ('+str(ClassMatchAll)+'%).\nThere were ' + str(SomeDiff) + 'different variant class predictions ('+str(ClassMismatchSome)+'%).\n')
print('VEP Impact\tSnpEff Impact\tNumber of Occurences',file=output)
for key in differences.keys():
	first,second = key.split('_')
	occurence = differences[key]
	print(first,second,occurence,sep='\t',file=output)



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

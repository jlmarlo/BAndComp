###recreating the script to maybe clean it up somewhat
import argparse
import json

###create arguments
parser = argparse.ArgumentParser(description='Comparing annotations between VEP, SnpeEff, and ANNOVAR')
parser.add_argument('annotation_file', type=str, help='Path to file containing a table of annotations from VEP, SNPEff, and ANNOVAR')
parser.add_argument('binning_dict', type=str, help='Path to json file containing common names (broad categories) and all variant categories associated with each common name')

parser.add_argument('vep_precedence', type=str, help='Path to json file indicating precedence of VEP variant classes')

parser.add_argument('eff_precedence',type=str,help='Path to json file indicating precedence of SnpEff variant classes')

args=parser.parse_args()


#####Read in input files
table = open(args.annotation_file, 'rt')
with open(args.binning_dict,'rt') as dictFile:
    classdict = json.load(dictFile)
with open(args.vep_precedence, 'rt') as vepFile:
    VepPre = json.load(vepFile)
with open(args.eff_precedence, 'rt') as effFile:
    EffPre = json.load(effFile)


##Create Output Files
final = open('outputs/final/FinalReport.txt','wt')

mismatchFile = open('outputs/final/Mismatches.txt','wt') 
V_Smismatch = open('outputs/final/V_SMismatch.txt','wt')
V_Amismatch = open('outputs/final/V_AMismatch.txt','wt')
S_Amismatch = open('outputs/final/S_AMismatch.txt','wt')

impactmismatches = open('outputs/final/ImpactMismatches.txt','wt')
classificationmismatches = open('outputs/final/ClassificationMismatches.txt','wt')
#fileName1 = open('outputs/final/VEPUniqueTranscripts.txt','wt')
#fileName2 = open('outputs/final/SnpEffUniqueTranscripts.txt','wt')

###Files for testing purposes
classes = open('outputs/final/ClassesList.txt','wt')
impacts = open('outputs/final/ImpactsList.txt','wt')
test = open('outputs/final/TestOutput.txt','wt')
special = open('outputs/final/OtherOutputs.txt','wt')

uneven = open('outputs/final/UniqueTranscriptsStats.txt','wt')

######Tracking Variable Creation##############
linenumber =1 #keep track of line of annotation table
astericks = 0 
nonastericks = 0
EffmultinoVEP = 0
VEPmultinoEff = 0
VEPCount = 0
EffCount = 0
impactMatch = 0
impactMismatch = 0
vepclasses = {}
effclasses = {}
annclasses = {}


VeptransC = 0
EfftransC = 0

table.readline() #remove header of table

##Name Change Code
def change_names(specName,tool):
    namechangecount = 0
    if isinstance(specName, list):
        commonFinal=[]
        for item in specName:
            MatchFound = False
            for common,namePos in classdict.items():
                if item in namePos:
                    commonFinal.append(common)
                    MatchFound = True
                    namechangecount +=1
                    break
            if MatchFound == False:
                print(f"No match found for {tool} class {item}")
    else:
        for common,namePos in classdict.items():
                if specName in namePos:
                    commonFinal=common
                    MatchFound = True
                    namechangecount+=1
                    break
        if MatchFound == False:
            print(f"No match found for {tool} class {item}")
    return commonFinal,namechangecount




def determine_precedence(AnnotationDictionary,tool):
    highestprecedence =100
    precidentClass = 'None'
    if tool =='VEP':
        precedenceList = VepPre 
    if tool =='SnpEff':
        precedenceList = EffPre
    for transcript in AnnotationDictionary:
        classific = AnnotationDictionary[transcript]['class']
        if isinstance(classific,list):
            for item in classific:
                if item in precedenceList:
                    newprecedence = precedenceList[item]
                    if newprecedence < highestprecedence:
                        highestprecedence = newprecedence
                        precidentClass = item
                else:
                    print(f'There is no precidence for {item}',file = test)
            AnnotationDictionary[transcript]['precedence']=highestprecedence
            AnnotationDictionary[transcript]['precClass']=precidentClass
            AnnotationDictionary[transcript]['precClassCommon'],fCount=change_names(precidentClass,tool)
        else:
            if classific in precedenceList:
                precedence = precedenceList[classific]
                AnnotationDictionary[transcript]['precedence'] = precedence
                AnnotationDictionary[transcript]['precClass'] = classific
                AnnotationDictionary[transcript]['precClassCommon'],fCount=change_names(classific,tool)
            else:
                print(f'There is no precidence for {classific}',file=test)
    return AnnotationDictionary



exceptions = ['5_prime_UTR_truncation&exon_loss_variant', '3_prime_UTR_truncation&exon_loss_variant']

def make_dictionary(annotation,tool):
    ## making dictionaries for VEP and snpEff. ANNOVAR is processed separately
    errorCount=0
    AnnDict = {}
    dupNumber = 1
    className = 'class'
    impactName = 'impact'
    number = 0
    prints='No'
    changecount = 0
    if ',' in annotation:
        multiAnno = annotation.split(',')
        for version in multiAnno:
            number +=1
            pieces = version.split('|')
            if pieces[6] == '':
                transcript = 'NoTranscript'
            elif '.' in pieces[6]:
                transcript = pieces[6].split('.')[0]
            else:
                transcript = pieces[6]
            if transcript in AnnDict:
                dupNumber +=1
                transcript = transcript +'^'+str(dupNumber)
                prints= 'Yes'
            AnnDict[transcript]={className:pieces[1],impactName:pieces[2]}
            if '&' in AnnDict[transcript][className]:
                number +=AnnDict[transcript][className].count('&')
                classes = AnnDict[transcript][className].split('&')
                merged_classes=[]
                i=0
                while i < len(classes): ##there are two cases where & is in classification.
                    ##This section recombines those classes together for faithful representation.
                    if i < len(classes)-1 and f"{classes[i]}&{classes[i+1]}" in exceptions:
                        merged_classes.append(f"{classes[i]}&{classes[i+1]}")
                        i+=2
                        number = number -1
                    else:
                        merged_classes.append(classes[i])
                        i+=1
                AnnDict[transcript][className] = merged_classes
            AnnDict[transcript]['common'],change=change_names(AnnDict[transcript][className],tool)
            changecount += change
    else:
        number = 1
        pieces = annotation.split('|')
        if pieces[6] == '':
            transcript = 'NoTranscript'
        elif '.' in pieces[6]:
            transcript = pieces[6].split('.')[0]
        else:
            transcript = pieces[6]
        if transcript in AnnDict:
                dupNumber +=2
                transcript = transcript +'^'+str(dupNumber)
                prints='Yes'
        AnnDict[transcript]={className:pieces[1],impactName:pieces[2]}
        if '&' in AnnDict[transcript][className]:
            number += AnnDict[transcript][className].count('&')
            classes = AnnDict[transcript][className].split('&')
            merged_classes=[]
            i=0
            while i < len(classes):
                if i < len(classes)-1 and f"{classes[i]}&{classes[i+1]}" in exceptions:
                    merged_classes.append(f"{classes[i]}&{classes[i+1]}")
                    i+=2
                    number = number -1
                else:
                    merged_classes.append(classes[i])
                    i+=1
            AnnDict[transcript][className] = merged_classes
        AnnDict[transcript]['common'],change = change_names(AnnDict[transcript][className],tool)
        changecount +=change
    AnnDictFinal = determine_precedence(AnnDict,tool)
    if changecount!=number:
        print(f"The number of name changes is {changecount} but the number of annotations is {number}\n",AnnDictFinal,file=test)
    return AnnDictFinal, number, changecount


def compareClass(annotation1, annotation2,impact,tool1,tool2):
    classMatch = 0
    classMismatch = 0
    impactMatch = 0
    impactMismatch = 0
    Ann1uni=0
    Ann2uni=0
    sharedtrans=0
    partialmatches =[]
    partialmatchtypes = []
    compCount =0
    partCount1 = 0
    partCount2 = 0
    imperrorType = 'none'
    classError = 'none'
    if len(annotation1) ==1 and len(annotation2) == 1:
        key1,value1 = next(iter(annotation1.items()))
        key2,value2 = next(iter(annotation2.items()))
        class1 = annotation1[key1]['common']
        class2 = annotation2[key2]['common']
        if isinstance(class1,list):
            partCount1 +=len(class1)
        else:
            partCount1 +=1
        if isinstance(class2, list):
            partCount2 +=len(class2)
        else:
            partCount2 +=1
        compCount =1
        if class1 == class2:
            classMatch +=1
        elif isinstance(class1,list) or isinstance(class2,list):
            if isinstance (class1, str):
                class1 = [class1]
            elif isinstance (class2, str):
                class2 = [class2]
            overlap = list(set(class1) & set(class2))
            class1Only = [x for x in class1 if x not in class2]
            class2Only = [x for x in class2 if x not in class1]

            unique = len(class1Only)+len(class2Only)
            pMatch = 100*(len(overlap)/(unique+len(overlap)))
            if pMatch == 0:
                classMismatch +=1
            elif pMatch ==100:
                classMatch +=1
                #print('This is wrong')
                #print(f'{class1}\t{class2}')
            else:
                PartialMatchType = str(len(class1Only))+'-'+str(len(class2Only))
                partialmatches.append(pMatch)
                partialmatchtypes.append(PartialMatchType)
        else:
            classMismatch +=1
            classError = class1+'-'+class2
            #print(f"Annotation1 classification: {annotation1[key1]['common']}. Annotation2 classification: {annotation2[key2]['common']}",file=classificationmismatches)
        if impact == 'Yes':
            if annotation1[key1]['impact'] == annotation2[key2]['impact']:
                impactMatch += 1
            else:
                impactMismatch +=1
                imperrorType = annotation1[key1]['impact'] + '-' +annotation2[key2]['impact']
                #print(imperrorType,file=test)
                #print(f"Annotation1 impact: {annotation1[key1]['impact']} Annotation2 impact: {annotation2[key2]['impact']}",file=impactmismatches)
        sharedtrans = 1
        if impact == 'Yes':
            return classMatch, classMismatch, impactMatch, impactMismatch, Ann1uni,Ann2uni,sharedtrans,partialmatches,partialmatchtypes,compCount,partCount1,partCount2,imperrorType,classError
        else:
            return classMatch, classMismatch, Ann1uni,Ann2uni,sharedtrans,partialmatches,partialmatchtypes,compCount,partCount1,partCount2,classError
    
    else:
        Unique1Classes = []
        Unique2Classes = []
        Ann1unique = set(annotation1.keys()) - set(annotation2.keys())
        Ann2unique = set(annotation2.keys()) - set(annotation1.keys())
        classError=[]
        if tool2 =='Annovar':
            ##iterate through all of the transcripts that the other annotation has
            for key in annotation1.keys():
            ##compare annovar class to annotatoin2 class
                class1 = annotation1[key]['common']
                key2,value2 = next(iter(annotation2.items()))
                class2 = annotation2[key2]['common']
                if class1 == class2:
                    classMatch +=1
                elif isinstance(class1,list) or isinstance(class2,list):
                    if isinstance (class1, str):
                        class1 = [class1]
                    elif isinstance (class2, str):
                        class2 = [class2]
                    overlap = list(set(class1) & set(class2))
                    class1Only = [x for x in class1 if x not in class2]
                    class2Only = [x for x in class2 if x not in class1]

                    unique = len(class1Only)+len(class2Only)
                    pMatch = 100*(len(overlap)/(unique+len(overlap)))
                    if pMatch == 0:
                        classMismatch +=1
                        classError.append('list_mismatch')
                    elif pMatch ==100:
                        classMatch +=1
                    else:
                        partialMatchType = str(len(class1Only))+'-'+str(len(class2Only))
                        partialmatches.append(pMatch)
                        partialmatchtypes.append(partialMatchType)
                else:
                    classMismatch +=1
                    classError.append(class1+'-'+class2)
        if tool1 == 'VEP' and tool2 == 'SnpEff':
            if len(Ann1unique) >0:
                for unique in Ann1unique:
                    #print(f"{unique} : {annotation1[unique]['class']} : {annotation1[unique]['common']}",file=fileName1)
                    Unique1Classes.append(unique+':'+str(annotation1[unique]['common']))
            if len(Ann2unique) >0:
                for unique in Ann2unique:
                    #print(f"{unique} : {annotation2[unique]['class']} : {annotation2[unique]['common']}",file = fileName2)
                    Unique2Classes.append(unique+':'+str(annotation2[unique]['common']))
        Ann1uni +=len(Ann1unique)
        Ann2uni +=len(Ann2unique)
        for item in Ann1unique:
            if isinstance(annotation1[item]['common'],list):
                partCount1+= len(annotation1[item]['common'])
            else:
                partCount1+=1
        for item in Ann2unique:
            if isinstance(annotation2[item]['common'],list):
                partCount2+=len(annotation2[item]['common'])
            else:
                partCount2+=1
        Shared = set(annotation1.keys()) & set(annotation2.keys())
        sharedtrans += len(Shared)
        imperrorType = []
        for transcript in Shared:
            compCount +=1
            class1 = annotation1[transcript]['common']
            class2 = annotation2[transcript]['common']
            if isinstance(class1,list):
                partCount1 +=len(class1)
            else:
                partCount1 +=1
            if isinstance(class2, list):
                partCount2 +=len(class2)
            else:
                partCount2 +=1
            if class1 == class2:
                classMatch +=1
            elif isinstance(class1,list) or isinstance(class2,list):
                if isinstance (class1, str):
                    class1 = [class1]
                elif isinstance (class2, str):
                    class2 = [class2]
                overlap = list(set(class1) & set(class2))
                class1Only = [x for x in class1 if x not in class2]
                class2Only = [x for x in class2 if x not in class1]

                unique = len(class1Only)+len(class2Only)
                pMatch = 100*(len(overlap)/(unique+len(overlap)))
                if pMatch == 0:
                    classMismatch +=1
                    classError.append('list_mismatch')
                elif pMatch ==100:
                    classMatch +=1
                else:
                    partialMatchType = str(len(class1Only))+'-'+str(len(class2Only))
                    partialmatches.append(pMatch)
                    partialmatchtypes.append(partialMatchType)
            else:
                classMismatch +=1
                classError.append(class1+'-'+class2)
            if impact == 'Yes':
                if annotation1[transcript]['impact'] == annotation2[transcript]['impact']:
                    impactMatch += 1
                    imperrorType.append('none')
                else:
                    impactMismatch +=1
                    imperrorType.append(annotation1[transcript]['impact']+'-'+annotation2[transcript]['impact'])
                    #print(imperrorType,file=test)
        if impact == 'Yes':
            return classMatch, classMismatch, impactMatch, impactMismatch, Ann1uni,Ann2uni,sharedtrans,partialmatches,partialmatchtypes,compCount,partCount1,partCount2,imperrorType,classError,Unique1Classes,Unique2Classes
        else:
            return classMatch, classMismatch, Ann1uni,Ann2uni,sharedtrans, partialmatches,partialmatchtypes,compCount,partCount1,partCount2,classError,Unique1Classes,Unique2Classes

def precedentCompare(annotation1, annotation2,Annovar):
    #find hihgest precident in dictionarys
    mismatch = 'none'
    impmatch = 'none'
    if Annovar == 'No':
        highestprecident = 100
        highestprecident2 = 100
        for item in annotation1:
            if annotation1[item]['precedence'] < highestprecident:
                precidentClass1 = annotation1[item]['precClassCommon']
                highestprecident = annotation1[item]['precedence']
                precimp1 = annotation1[item]['impact']
        for item in annotation2:
            if annotation2[item]['precedence'] < highestprecident2:
                precidentClass2 = annotation2[item]['precClassCommon']
                highestprecident2 = annotation2[item]['precedence']
                precimp2 = annotation2[item]['impact']
        if precidentClass1 == precidentClass2:
            match = 'Yes'
        else:
            match = 'No'
            #print(precidentClass1, precidentClass2,file=test)
            mismatch = precidentClass1+'-'+precidentClass2
        if precimp1 == precimp2:
            impmatch = 'Yes'
        else:
            impmatch = 'No'
    else:
        highestprecident = 100
        for item in annotation1:
            if annotation1[item]['precedence'] < highestprecident:
                precidentClass = annotation1[item]['precClassCommon']
                highestprecident = annotation1[item]['precedence']
        key1,value1 = next(iter(annotation2.items()))
        precidentClass2 = annotation2[key1]['common']
        if isinstance(precidentClass2, list):
            if precidentClass in precidentClass2:
                match = 'Yes'
            else:
                match = 'No'
        else:
            if precidentClass == precidentClass2:
                match ='Yes'
            else:
                match ='No'
                mismatch = precidentClass+'-'+precidentClass2
    return match,mismatch,impmatch

def track_changes(variant_dict):
    changes = {}
    impacts = set()
    classifications= set()
    multiTransc = 'No'
    if len(variant_dict) > 1:
        multiTransc ='Yes'
    for transcript in variant_dict.keys():
        impacts.add(variant_dict[transcript]['impact'])
        classifications.add(str(variant_dict[transcript]['class']))

        # Initialize a dictionary to store the results for each variant
    changes = {
            'impact': {
                'unique_count': len(impacts),
                'unique_values': list(impacts),
                'changes_count': len(impacts) - 1 if len(impacts) > 1 else 0,
            },
            'classification': {
                'unique_count': len(classifications),
                'unique_values': list(classifications),
                'changes_count': len(classifications) - 1 if len(classifications) > 1 else 0,
            }
        }

    return changes,multiTransc

def count_changes(variant_results):
    impact_change_count = 0
    impact_multichange_count = 0
    classification_change_count = 0
    impact_nochange = 0
    classification_nochange = 0
    classification_multichange_count=0
    for variant, details in variant_results.items():
        # Check if there's a change in impact
        if details['impact']['changes_count'] > 1:
            impact_multichange_count += 1
        elif details['impact']['changes_count'] >0:
            impact_change_count+=1
        else:
            impact_nochange +=1

        # Check if there's a change in classification
        if details['classification']['changes_count'] > 1:
            classification_multichange_count += 1
        elif details['classification']['changes_count'] >0:
            classification_change_count+=1
        else:
            classification_nochange +=1
    return impact_change_count, classification_change_count, impact_nochange, classification_nochange, classification_multichange_count, impact_multichange_count


VEPChangesDictionary = {}
EffChangesDictionary = {}


V_SMatch = 0
V_SMismatch = 0
V_SPartial = 0
V_SComparisons = 0
impactMatch = 0
impactMismatch = 0
VepPartCount = 0
V_AMatch = 0
V_AMismatch = 0
V_APartial = 0
V_AComparisons = 0
S_AMatch = 0
S_AMismatch = 0
S_APartial = 0
S_AComparisons = 0


sharedtrans=0
VEPtransc = 0
Efftransc = 0
TotalVepNum = 0
TotalEffNum=0
VepErrors = 0
EffErrors = 0
Effchanges = 0
Annchanges = 0
VEPchanges = 0
impactMatch=0
impactMismatch =0
AnnDictionary={}
singleannotation=0
functionsingle = 0
NumberComparisons = 0
VepPartCount = 0
EffPartCount = 0
V_SprecidentMatch = 0
V_SprecidentMismatch = 0
V_AprecidentMatch = 0
V_AprecidentMismatch = 0
S_AprecidentMatch = 0
S_AprecidentMismatch = 0
V_SPrecDiff = {}
V_APrecDiff={}
S_APrecDiff={}
V_SprecimpMatch = 0
V_SprecimpMismatch = 0

effimpacts = {}
vepimpacts = {}
impactDiffDict ={}
V_SclassDiffDict = {}
V_AclassDiffDict = {}
S_AclassDiffDict = {}
VUniqueDict={}
EUniqueDict={}
V_SPartialDict = {}
V_APartialDict = {}
S_APartialDict = {}
VepUniqueNoTranscript = 0
VepActualUniTranscript = 0
BetweenGenes = 0
duplicateTranscript = 0
EffActualUniTranscript = 0
EffGene = 0
VepChangesDict = {}
EffChangesDict = {}
V_SPartsDict = {}
V_APartsDict = {}
S_APartsDict = {}
VEPmultiCount = 0
VEPsingleCount = 0
EffmultiCount = 0
EffsingleCount = 0
for variant in table:
    linenumber +=1
    variant = variant.strip()
    infoList = variant.split('\t')
    if '*' in infoList[3]:
        astericks +=1
        continue
    variantID = infoList[0] + '['+infoList[1]+']'+infoList[2]+'-'+infoList[3]
    VEPAnnotation = infoList[4] #save VEP annotation
    EffAnnotation = infoList[5] #save SnpEff annotation
    AnnovarAnnotation = infoList[6:11]
    
    #####Processing ANNOVAR annotations
    ann_area = AnnovarAnnotation[0]

    if '3b' in ann_area:
        ann_class = ann_area.split('\\x3b')
        if 'exonic' in ann_class:
            ann_class.append(AnnovarAnnotation[3])
    elif ann_area == 'exonic':
        ann_class = AnnovarAnnotation[3]
    else:
        ann_class = ann_area
    AnnDictionary['annovar']={'class':ann_class,'common':'none'}


    #######Making list of ann ANNOVAR classes#######
    if isinstance(ann_class,list):
        for item in ann_class:
            if item in annclasses:
                annclasses[item] +=1
            else:
                annclasses[item] =1

    else:
        if ann_class in annclasses:
            annclasses[ann_class]+=1
        else:
            annclasses[ann_class]=1

    ######Making Dictionary for ANNOVAR######
    AnnDictionary['annovar']['common'],namecount=change_names(ann_class,'Annovar')
    Annchanges +=namecount

    
    
    
    ######making lists for snpeff and vep classes

    VepDictionary,VepNumAnn,Vnamecount = make_dictionary(VEPAnnotation,"VEP")

    EffDictionary,EffNumAnn,Enamecount = make_dictionary(EffAnnotation,"SnpEff")
    

    ###Tracking when transcripts change classification and impact predictions
    VepChangesDict[variantID],VepmultiTranscr=track_changes(VepDictionary)
    EffChangesDict[variantID],EffmultiTranscr=track_changes(EffDictionary)

    if VepmultiTranscr == 'Yes':
        VEPmultiCount +=1
    else:
        VEPsingleCount +=1
    if EffmultiTranscr == 'Yes':
        EffmultiCount +=1
    else:
        EffsingleCount +=1
    ##trouble shooting variables###
    TotalVepNum += VepNumAnn
    TotalEffNum+= EffNumAnn
    Effchanges += Enamecount
    VEPchanges += Vnamecount


####Code to loop through the dictionaries and count up instances of each variant classification type
    for Vtranscript in VepDictionary.keys():
        classif = VepDictionary[Vtranscript]['class']
        impact = VepDictionary[Vtranscript]['impact']
        if isinstance(classif,list):
            for item in classif:
                if item in vepclasses:
                    vepclasses[item]+=1
                else:
                    vepclasses[item] =1
        else:
            if classif in vepclasses:
                vepclasses[classif]+=1
            else:
                vepclasses[classif] =1
        if impact in vepimpacts:
            vepimpacts[impact] +=1
        else:
            vepimpacts[impact]=1


    for Etranscript in EffDictionary.keys():
        classif = EffDictionary[Etranscript]['class']
        impact = EffDictionary[Etranscript]['impact']
        if isinstance(classif,list):
            for item in classif:
                if item in effclasses:
                    effclasses[item]+=1
                else:
                    effclasses[item] =1
        else:
            if classif in effclasses:
                effclasses[classif]+=1
            else:
                effclasses[classif]=1
        if impact in effimpacts:
            effimpacts[impact]+=1
        else:
            effimpacts[impact]=1

    ######SnpEff and VEP transcript comparisons######

    V_SCompStats = compareClass(VepDictionary,EffDictionary,'Yes','VEP','SnpEff')

    V_ACompStats = compareClass(VepDictionary,AnnDictionary,'No','VEP','Annovar')

    S_ACompStats = compareClass(EffDictionary,AnnDictionary,'No','SnpEff','Annovar')


    ####Counting Up Stats

    V_SMatch +=V_SCompStats[0]
    V_SMismatch +=V_SCompStats[1]
    V_SPartial +=len(V_SCompStats[7])
    V_SPart = V_SCompStats[7]
    V_SPartType = V_SCompStats[8]
    V_SComparisons +=V_SCompStats[9]
    impactMatch += V_SCompStats[2]
    impactMismatch +=V_SCompStats[3]
    VepPartCount += V_SCompStats[10]
    EffPartCount += V_SCompStats[11]
    impactDifferences = V_SCompStats[12]
    V_SclassDifferences = V_SCompStats[13]
    
    #####disastrously analysing unique transcript contents
    try:
        VUnique = V_SCompStats[14]
        if len(VUnique) > 0:
            for item in VUnique:
                transcript,classtype = item.split(':')
                if transcript == 'NoTranscript':
                    VepUniqueNoTranscript +=1
                    #print(variant,file=test)
                elif 'T' in transcript:
                    VepActualUniTranscript +=1
                    #print(variant,file=test)
                else:
                    print(f'VEP unique transcript: {transcript}')
                if ',' in classtype:
                    classtype = classtype.rstrip("']").lstrip("['")
                    pieces = classtype.split(',')
                    for piece in pieces:
                        piece = piece.rstrip("'")
                        piece = piece.lstrip()
                        piece = piece.lstrip("'")
                        if piece in VUniqueDict:
                            VUniqueDict[piece]+=1
                        else:
                            VUniqueDict[piece]=1
                else:
                    if classtype in VUniqueDict:
                        VUniqueDict[classtype]+=1
                    else:
                        VUniqueDict[classtype]=1
        EUnique = V_SCompStats[15]
        if len(EUnique) > 0:
            for item in EUnique:
                transcript,classtype = item.split(':')
                if '-' in transcript:
                    BetweenGenes +=1
                    #print(classtype,file=test)
                elif '^' in transcript:
                    duplicateTranscript +=1
                    #print(classtype,file=test)
                elif 'T' in transcript:
                    EffActualUniTranscript +=1
                elif 'G' in transcript:
                    EffGene +=1
                    print(classtype,file=test)
                else:
                    print(f'Eff unique transcript: {transcript}')
                if ',' in classtype:
                    classtype = classtype.rstrip("']").lstrip("['")
                    pieces = classtype.split(',')
                    for piece in pieces:
                        piece = piece.rstrip("'")
                        piece = piece.lstrip()
                        piece = piece.lstrip("'")
                        if piece in EUniqueDict:
                            EUniqueDict[piece]+=1
                        else:
                            EUniqueDict[piece]=1
                else:
                    if classtype in EUniqueDict:
                        EUniqueDict[classtype]+=1
                    else:
                        EUniqueDict[classtype]=1
    except IndexError:
        pass


    V_AMatch +=V_ACompStats[0]
    V_AMismatch +=V_ACompStats[1]
    V_APartial +=len(V_ACompStats[5])
    V_AComparisons +=V_ACompStats[7]
    V_AclassDifferences = V_ACompStats[10]
    V_APart = V_ACompStats[5]
    V_APartType = V_ACompStats[6]

    S_AMatch +=S_ACompStats[0]
    S_AMismatch +=S_ACompStats[1]
    S_APartial +=len(S_ACompStats[5])
    S_APart = S_ACompStats[5]
    S_AComparisons +=S_ACompStats[7]
    S_AclassDifferences = S_ACompStats[10]
    S_APartType = S_ACompStats[6]

    ######counting up partial match percentages#######
    for percent in V_SPart:
        if percent in V_SPartialDict:
            V_SPartialDict[percent]+=1
        else:
            V_SPartialDict[percent]=1
    
    for kind in V_SPartType:
        if kind in V_SPartsDict:
            V_SPartsDict[kind]+=1
        else:
            V_SPartsDict[kind]=1

    for percent in V_APart:
        if percent in V_APartialDict:
            V_APartialDict[percent]+=1
        else:
            V_APartialDict[percent]=1

    for kind in V_APartType:
        if kind in V_APartsDict:
            V_APartsDict[kind]+=1
        else:
            V_APartsDict[kind]=1

    for percent in S_APart:
        if percent in S_APartialDict:
            S_APartialDict[percent]+=1
        else:
            S_APartialDict[percent]=1

    for kind in S_APartType:
        if kind in S_APartsDict:
            S_APartsDict[kind]+=1
        else:
            S_APartsDict[kind]=1



    ####counting up number of impact disagreements
    if isinstance(impactDifferences,list):
        for item in impactDifferences:
            if item in impactDiffDict:
                impactDiffDict[item] +=1
            else:
                impactDiffDict[item]=1
    else:
        if impactDifferences in impactDiffDict:
            impactDiffDict[impactDifferences] +=1
        else:
            impactDiffDict[impactDifferences]=1
    


    ####counting up number of vep snpeff classification disagreements
    if isinstance(V_SclassDifferences,list):
        for thing in V_SclassDifferences:
            if thing in V_SclassDiffDict:
                V_SclassDiffDict[thing]+=1
            else:
                V_SclassDiffDict[thing]=1
    else:
        if V_SclassDifferences in V_SclassDiffDict:
            V_SclassDiffDict[V_SclassDifferences]+=1
        else:
            V_SclassDiffDict[V_SclassDifferences]=1


    ####counting up number of vep annovar classification disagreements
    if isinstance(V_AclassDifferences,list):
        for thing in V_AclassDifferences:
            if thing in V_AclassDiffDict:
                V_AclassDiffDict[thing]+=1
            else:
                V_AclassDiffDict[thing]=1
    else:
        if V_AclassDifferences in V_AclassDiffDict:
            V_AclassDiffDict[V_AclassDifferences]+=1
        else:
            V_AclassDiffDict[V_AclassDifferences]=1


    ###counting up number of snpeff annovar classification disagreements
    if isinstance(S_AclassDifferences,list):
        for thing in S_AclassDifferences:
            if thing in S_AclassDiffDict:
                S_AclassDiffDict[thing]+=1
            else:
                S_AclassDiffDict[thing]=1
    else:
        if S_AclassDifferences in S_AclassDiffDict:
            S_AclassDiffDict[S_AclassDifferences]+=1
        else:
            S_AclassDiffDict[S_AclassDifferences]=1


    ####counting up matches when using precidences
    ###Vep - Snpeff
    matchStat,mismatch,precimpMatch=precedentCompare(VepDictionary,EffDictionary,'No')
    if matchStat == 'Yes':
        V_SprecidentMatch +=1
    else:
        V_SprecidentMismatch +=1
    
    if mismatch in V_SPrecDiff:
        V_SPrecDiff[mismatch]+=1
    else:
        V_SPrecDiff[mismatch]=1

    if precimpMatch == 'Yes':
        V_SprecimpMatch +=1
    else:
        V_SprecimpMismatch +=1

    ###Vep-Annovar
    matchStat,mismatch,precimpMatch=precedentCompare(VepDictionary,AnnDictionary,'Yes')
    if matchStat == 'Yes':
        V_AprecidentMatch +=1
    else:
        V_AprecidentMismatch +=1

    if mismatch in V_APrecDiff:
        V_APrecDiff[mismatch]+=1
    else:
        V_APrecDiff[mismatch]=1

    ##Eff - Annovar
    matchStat,mismatch,precimpMatch=precedentCompare(EffDictionary,AnnDictionary,'Yes')
    if matchStat == 'Yes':
        S_AprecidentMatch +=1
    else:
        S_AprecidentMismatch +=1

    if mismatch in S_APrecDiff:
        S_APrecDiff[mismatch]+=1
    else:
        S_APrecDiff[mismatch]=1

EffTranscriptChanges = count_changes(EffChangesDict)
VepTranscriptChanges = count_changes(VepChangesDict)

print(f"{VEPmultiCount} <- VEP multi. {VEPsingleCount} <- VEP single. {EffmultiCount} <- Eff multi. {EffsingleCount} <- Eff single.")

########################MAKING FINAL OUTPUTS###########################
print(f"{astericks} lines were removed due to invalid alternate alleles (*)\n\n",file=final)
print("PERCENT MATCHES OF TOOLS\n",file=final)
print("VEP versus SnpEff",file=final)
print(f"Matches: {V_SMatch} ({100*(V_SMatch/(V_SMismatch+V_SPartial+V_SMatch))}%). \nMismatches : {V_SMismatch} ({100*(V_SMismatch/(V_SMismatch+V_SPartial+V_SMatch))}%). \nPartial matches : {V_SPartial} ({100*(V_SPartial/(V_SMismatch+V_SPartial+V_SMatch))}%)",file=final)
print(f"Total: {V_SMatch+V_SMismatch+V_SPartial}",file=final)
print(f"Total Theoretical comparison count: {V_SComparisons}\n",file=final)
print(f"Impact matches : {impactMatch} ({100*(impactMatch/(impactMatch+impactMismatch))}%). \nTotal mismatches: {impactMismatch} ({100*(impactMismatch/(impactMatch+impactMismatch))}%).\nTotal = {impactMatch + impactMismatch}\n",file=final)

#print(f"Total number of parts in comparisons. VEP: {VepPartCount}. Eff: {EffPartCount}.")

print('VEP versus ANNOVAR',file=final)
print(f"Matches : {V_AMatch} ({100*(V_AMatch/(V_AMismatch+V_APartial+V_AMatch))}%). \nMismatches: {V_AMismatch} ({100*(V_AMismatch/(V_AMismatch+V_APartial+V_AMatch))}%). \nPartial : {V_APartial} ({100*(V_APartial/(V_AMismatch+V_APartial+V_AMatch))}%)",file=final)
print(f"Total: {V_AMatch+V_AMismatch+V_APartial}",file=final)
print(f"Total Theoretical comparison count: {V_AComparisons}\n",file=final)

print('SnpEff versus ANNOVAR',file=final)
print(f"Matches: {S_AMatch} ({100*(S_AMatch/(S_AMismatch+S_APartial+S_AMatch))}%). \nMismatches: {S_AMismatch} ({100*(S_AMismatch/(S_AMismatch+S_APartial+S_AMatch))}%). \nPartial : {S_APartial} ({100*(S_APartial/(S_AMismatch+S_APartial+S_AMatch))}%)",file=final)
print(f"Total: {S_AMatch+S_AMismatch+S_APartial}",file=final)
print(f"Total Theoretical comparison count: {S_AComparisons}\n\n",file=final)


print('---------------------PRECIDENT COMPARISONS---------------------\n',file=final)
print(f'V_S matches: {V_SprecidentMatch} ({100*(V_SprecidentMatch/(V_SprecidentMatch+V_SprecidentMismatch))}%). Mismatches: {V_SprecidentMismatch} ({100*(V_SprecidentMismatch/(V_SprecidentMatch+V_SprecidentMismatch))}%). Total = {V_SprecidentMatch+V_SprecidentMismatch}',file=final)
print(f'V_A matches: {V_AprecidentMatch} ({100*(V_AprecidentMatch/(V_AprecidentMatch+V_AprecidentMismatch))}%). Number of mismatches {V_AprecidentMismatch} ({100*(V_AprecidentMismatch/(V_AprecidentMatch+V_AprecidentMismatch))}%).  Total = {V_AprecidentMatch+V_AprecidentMismatch}',file=final)
print(f'S_A matches: {S_AprecidentMatch} ({100*(S_AprecidentMatch/(S_AprecidentMatch+S_AprecidentMismatch))}%). Number of mismatches {S_AprecidentMismatch} ({100*(S_AprecidentMismatch/(S_AprecidentMatch+S_AprecidentMismatch))}%). Total = {S_AprecidentMatch+S_AprecidentMismatch}',file=final)


print(f'Impact matches: {V_SprecimpMatch} ({100*(V_SprecimpMatch/(V_SprecimpMatch+V_SprecimpMismatch))}%). Number of mismatches {V_SprecimpMismatch} ({100*(V_SprecimpMismatch/(V_SprecimpMatch+V_SprecimpMismatch))}%). Total = {V_SprecimpMatch+V_SprecimpMismatch}',file=final)
####Closing Files

Vtotal = 0
Etotal =0
ATotal = 0
#########Creating Outputs of Class Lists#######
print("ANNOVAR CLASSES",file=classes)
for key,value in annclasses.items():
    print(f"{key} : {value}",file = classes)
    ATotal += value
print("\n\nVEP CLASSES", file = classes)
for key,value in vepclasses.items():
    print(f"{key} : {value}",file = classes)
    Vtotal += value
print("\n\nEFF CLASSES", file = classes)
for key,value in effclasses.items():
    print(f"{key} : {value}",file = classes)
    Etotal += value

print(f'\nAnnovar total classes: {ATotal}. VEP total classes: {Vtotal}. SnpEff total classes: {Etotal}.',file = classes)

print('VEP Impacts',file=impacts)
for key,value in vepimpacts.items():
    print(f'{key} : {value}',file=impacts)
print('\n\nSnpEff Impacts',file = impacts)
for key,value in effimpacts.items():
    print(f'{key} : {value}',file=impacts)


print('Impact Mismatches\n',file=mismatchFile)
for key,value in impactDiffDict.items():
    print(f'{key}\t {value}',file=mismatchFile)
print('\n\nVEP\tSnpEff',file=mismatchFile)
for key,value in V_SclassDiffDict.items():
    print(f'{key}\t{value}',file=mismatchFile)
print('\n\nVEP\tAnnovar',file=mismatchFile)
for key,value in V_AclassDiffDict.items():
    print(f'{key}\t{value}',file=mismatchFile)
print('\n\nSnpEff\tAnnovar',file=mismatchFile)
for key,value in S_AclassDiffDict.items():
    print(f'{key}\t{value}',file=mismatchFile)

print('-------------------Precident Differences-------------------',file=mismatchFile)

print('VEP SNPEff',file=mismatchFile)
for key,value in V_SPrecDiff.items():
    print(f'{key}\t{value}',file=mismatchFile)
print('\n\nVEP ANNOVAR',file=mismatchFile)
for key,value in V_APrecDiff.items():
    print(f'{key}\t{value}',file=mismatchFile)
print('\n\nSnpEff ANNOVAR',file=mismatchFile)
for key,value in S_APrecDiff.items():
    print(f'{key}\t{value}',file=mismatchFile)


print('PARTIAL VALUES DICTIONARY TEST',file=final)
print('\nVEP SnpEff Partials',file=final)
for key,value in V_SPartialDict.items():
    print(f'{key}\t{value}',file=final)
print('\nVEP ANNOVAR partials',file=final)
for key,value in V_APartialDict.items():
    print(f'{key}\t{value}',file=final)
print('\nSnpEff ANNOVAR partials',file=final)
for key,value in S_APartialDict.items():
    print(f'{key}\t{value}',file=final)


print('PARTIAL MATCHES TYPES DICTIONARY TEST',file=final)
print('\nVEP SnpEff Partials',file=final)
for key,value in V_SPartsDict.items():
    print(f'{key}\t{value}',file=final)
print('\nVEP ANNOVAR partials',file=final)
for key,value in V_APartsDict.items():
    print(f'{key}\t{value}',file=final)
print('\nSnpEff ANNOVAR partials',file=final)
for key,value in S_APartsDict.items():
    print(f'{key}\t{value}',file=final)







print('UNIQUE TRANSCRIPT CLASSES DICTIONARY TESTS\n\n',file=uneven)
print('VEP Dictionary',file=uneven)
for key,value in VUniqueDict.items():
    print(f'{key}\t{value}',file=uneven)
print("\nSnpEff Dictionary",file=uneven)
for key,value in EUniqueDict.items():
    print(f'{key}\t{value}',file=uneven)

print(f'\n\nVEP No Transcripts: {VepUniqueNoTranscript} \nVEP Unique Transcripts: {VepActualUniTranscript}. \n EFF Between Genes (Gene-Gene): Eff Unique Counts: {BetweenGenes} \nEff duplicate transcript: {duplicateTranscript} \nEff Unique Transcript: {EffActualUniTranscript} Eff Gene Transcript: {EffGene} \n\n',file=uneven)


##Transcript Changing Stats

print('----------------VEP transcript stats--------------',file=final)
print(f'Times impact prediction changed: {VepTranscriptChanges[0]+VepTranscriptChanges[5]}',file=final)
print(f'Times impact prediction changed more than once: {VepTranscriptChanges[5]}',file=final)
print(f'Time impact prediction did not change: {VepTranscriptChanges[2]}\n\n',file=final)
print(f'Time classificatoin changed: {VepTranscriptChanges[1]+VepTranscriptChanges[4]}',file=final)
print(f'Time classification changed more than once: {VepTranscriptChanges[4]}',file=final)
print(f'Time classification did not change: {VepTranscriptChanges[3]}\n\n',file=final)


print('-----------------SnpEff transcript stats-------------', file=final)
print(f'Times impact prediction changed: {EffTranscriptChanges[0]+EffTranscriptChanges[5]}',file=final)
print(f'Times impact prediction changed more than once: {EffTranscriptChanges[5]}',file=final)
print(f'Time impact prediction did not change: {EffTranscriptChanges[2]}\n\n',file=final)

print(f'Time classificatoin changed: {EffTranscriptChanges[1]+EffTranscriptChanges[4]}',file=final)
print(f'Time classification changed more than once: {EffTranscriptChanges[4]}',file=final)
print(f'Time classification did not change: {EffTranscriptChanges[3]}\n\n',file=final)

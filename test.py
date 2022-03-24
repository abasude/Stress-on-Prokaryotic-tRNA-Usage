genomeDict = {'YPK_0001': 'atcatc', 'YPK_0002': 'ctgatc', 'YPK_0005': 'atcctgctg'}

# iterate through genes
for gene, sequence in genomeDict.items():
    codonDict = {} # hold codon count for each gene
    totalCodonCount = 0  # keeps track of number of codons in each seq

    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3].upper() # create codon window
        

        # if the codon is not already in the dictionary, create key
        if codon not in codonDict:
            codonDict[codon] = 0

        codonDict[codon] += 1
        totalCodonCount += 1

    # print(totalCodonCount)
    # normalize codon values
    for codonKey, numCodons in codonDict.items():
        if totalCodonCount != 0:
            codonDict[codonKey] = numCodons / totalCodonCount

    genomeDict[gene] = codonDict

# print("genomeDict") # codonDictDict in real
# print(genomeDict)

import re

line = '>YPK_0001 [gene=YPK_0001] [protein=chromosomal replication initiator protein DnaA] [protein_id=ACA66315.1] [location=37..1389]'
header = re.split(" ", line[1:], 1)
header = header[0]
# print(header)

#,Geneid,Chr,Start,End,Strand,Length,sampleTreat_1.sam,sampleTreat_2.sam,sampleTreat_3.sam,sampleCtrl_1.sam,sampleCtrl_2.sam,sampleCtrl_3.sam,sampleTreat_1.sam.tpm,sampleTreat_2.sam.tpm,sampleTreat_3.sam.tpm,sampleCtrl_1.sam.tpm,sampleCtrl_2.sam.tpm,sampleCtrl_3.sam.tpm,sampleTreat_1.sam.cpm,sampleTreat_2.sam.cpm,sampleTreat_3.sam.cpm,sampleCtrl_1.sam.cpm,sampleCtrl_2.sam.cpm,sampleCtrl_3.sam.cpm
reads = ["0,YPK_0001,NC_010465,37,1389,+,1353,5,4,19,2,6,2,82.1296,64.0273,312.4351,40.1565,115.2275,48.8871,148.5134,121.4255,619.4777,76.7342,220.135,89.4815", 
         "1,YPK_0002,NC_010465,1394,2494,+,1101,7,12,8,8,14,6,141.2988,236.0462,161.6615,197.3905,330.4027,180.2297,207.9187,364.2766,260.8327,306.9368,513.6484,268.4444", 
         "2,YPK_0003,NC_010465,2667,3752,+,1086,0,7,0,0,1,0,0.0,139.5955,0.0,0.0,23.9262,0.0,0.0,212.4947,0.0,0.0,36.6892,0.0"]

expressionDict = {} 

for read in reads:
    splits = re.split(",", read)
    # [sampletpm, controltpm]
    expressionDict[splits[1]] = [splits[13], splits[16]]
   

# print("expressionDict")
# print(expressionDict)

weightedCodonDict = {}
for gene in genomeDict:
    if gene in expressionDict:
        codonDict = {}
        for codon, count in genomeDict[gene].items():
            codonDict[codon] = [float(expressionDict[gene][0]) * float(count),
                                float(expressionDict[gene][1]) * float(count)]

            weightedCodonDict[gene] = codonDict

# print("weightedCodonDict")
# print(weightedCodonDict)


import numpy as np
import matplotlib.pyplot as plt

# weightedCodonDict = {'YPK_0001': {'ATC': [82.1296, 40.1565]}, 'YPK_0002': {'CTG': [70.6494, 98.69525], 'ATC': [70.6494, 98.69525]}, 'YPK_0003': {'CTG': [70.6494, 98.69525], 'ATC': [70.6494, 98.69525]}}


codonDict = {} # holds added weighted codon scores across all genes 
for gene in weightedCodonDict:
    for codon, value in weightedCodonDict[gene].items():
        # add codon without gene to dictionary if not there
        if codon not in codonDict:
            codonDict[codon] = value
        # add the values for the codons if already there
        else:
            codonDict[codon] = [codonDict[codon][0] + value[0], codonDict[codon][1]+ value[1]] 


# print(codonDict)
# # plot bar graph 
# codons = list(codonDict.keys()) # list of codon names
# counts = list(codonDict.values()) # list of lists codon values for treatment and control
# treatment = [] # list of treatment weighted codon value
# control = [] # list of control weighted codon values

# # seperate treatment and control values
# for count in counts:
#     treatment.append(count[0])
#     control.append(count[1])

# # plot grouped bar 
# x = np.arange(len(codons))  # the label locations
# width = 0.2  # the width of the bars

# # fig, ax = plt.subplots()
# plt.bar(x - width/2, control, width, label='control', color = 'cyan' )
# plt.bar(x + width/2, treatment, width, label='treatment', color = 'red')

# plt.title('Weighted Codon Counts for Treatment VS Control Yersenia Pseudotuberculosis')
# plt.xticks(rotation='vertical')
# plt.xticks(x, codons)
# plt.xlabel("Codons")
# plt.ylabel("Scores")
# plt.legend(["Control","Treatment"])

# codonToAaDict = {
#     'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
#     'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
#     'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
#     'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
#     'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
#     'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
#     'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
#     'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
#     'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
#     'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
#     'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
#     'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
#     'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
#     'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
#     'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
#     'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

# aaColorDict = {
#     'M': 'gray', 'P': 'orange', 'A': 'green', 'V': 'red', 'L': 'purple',
#     'I': 'brown','G': 'blue', 'C': 'olive', 'F': 'pink', 'Y': 'cyan',		
#     'W': 'maroon',	'H': 'teal', 'K': 'indigo',	'R': 'lime', 'Q': 'magenta', 
#     'N': 'olive', 'E': 'springgreen', 'D': 'coral', 'S': 'mediumvioletred', 'T': 'gold'}

#     # color code the codons by aa
# for codonInd in range(len(codons)):
#     plt.gca().get_xticklabels()[codonInd].set_color(aaColorDict[codonToAaDict[codons[codonInd]]])

# plt.show()

# differences as keys and tpms as values
differenceTupList = []
for gene, tpmList in expressionDict.items():
    differenceTupList.append((gene, abs(float(tpmList[0]) - float(tpmList[1]))))

differenceTupList = sorted(differenceTupList, key=lambda x:x[1], reverse = True)

# print("gene\t\tTPM Difference")
# for i in range(0,2):
#     print("{}\t{}".format(differenceTupList[i][0], differenceTupList[i][1]))
codonDictDict = genomeDict



# find average codon use across all genes for each codon
# Ex. {'ATC': 0.611111111111111, 'CTG': 0.38888888888888884}
# {'codon': "average normalized codon use across all genes"}
averageCodonDict = {}
geneCount = 0 # keeps track of number of genes
for gene, codonDict in codonDictDict.items():
    geneCount += 1
    for codon in codonDict:
        # if codon not in the dict create codon key with normalized codon count value
        if codon not in averageCodonDict:
            averageCodonDict[codon] = codonDict[codon]
        # if in the dict add normalized codon count value
        else:
            averageCodonDict[codon] = averageCodonDict[codon] + codonDict[codon]

# divide by number of genes to get average codon count across all genes
for codon in averageCodonDict:
    averageCodonDict[codon] = averageCodonDict[codon] / geneCount



# print("averageCodonDict")
# print(averageCodonDict)


# find gene with highest disparity in average codon count for each codon across genes
# Ex. {'ATC': ['YPK_0001', 0.38888888888888895], 'CTG': ['YPK_0005', 0.2777777777777778]}
# {codon: [gene with highest difference from avg, difference]}
geneCodonDisparityDict = {}

# initialize with placeholder values
for codon in averageCodonDict:
    geneCodonDisparityDict[codon] = ['gene', 0]

for gene, codonDict in codonDictDict.items():
    for codon in codonDict: 
        # calculates difference between avg
        difference = abs(averageCodonDict[codon]-codonDict[codon])
        # if difference is bigger than what is already in dict
        if difference > geneCodonDisparityDict[codon][1]:
            # update dict
            geneCodonDisparityDict[codon] = [gene, difference]


geneCodonDisparityDict = sorted(geneCodonDisparityDict.items())
print("geneCodonDisparityDict")
print(geneCodonDisparityDict)

# print("codon\tgene\t\tdisparity")
# for codon in geneCodonDisparityDict:
#     print("{}\t{}\t{}".format(codon, geneCodonDisparityDict[codon][0], geneCodonDisparityDict[codon][1]))

# readStressed =["lcl|NZ_CP020549.1_cds_WP_000660615.1_1	1362	1112.000	40.101822	479.000", 
#          "lcl|NZ_CP020549.1_cds_WP_001208981.1_10	1278	1028.000	26.443754	292.000", 
#          "lcl|NZ_CP020549.1_cds_13	805	555.000	4.864502	29.000"]

# readControl = ["lcl|NZ_CP020549.1_cds_WP_000660615.1_1	1362	1112.000	40.101822	479.000", 
#          "lcl|NZ_CP020549.1_cds_WP_001208981.1_10	1278	1028.000	26.443754	292.000", 
#          "lcl|NZ_CP020549.1_cds_13	805	555.000	4.864502	29.000"]


# expressionDict = {} 

# for read in readStressed:
#     splits = re.split(" +|\t", read)
#     header = splits[0]
#     headerSplit =  re.split("_", header)
#     gene = "_".join(headerSplit[3:5])
#     expressionDict[gene] = splits[3]

# for read in readControl:
#     splits = re.split(" +|\t", read)
#     header = splits[0]
#     headerSplit =  re.split("_", header)
#     gene = "_".join(headerSplit[3:5])
#     if gene in expressionDict:
#         expressionDict[gene] = [expressionDict[gene], splits[3]]

# # special case because of weird output
# remove = []
# for gene in expressionDict:
#     if gene[0] != "W":
#         remove.append(gene)

# for gene in remove:
#     expressionDict.pop(gene)


# print("expressionDict")
# print(expressionDict)

# header =  ">lcl|NZ_CP020549.1_cds_WP_000660615.1_1 [gene=dnaA] [locus_tag=SPNHU17_RS00005] [db_xref=GeneID:66805161] [protein=chromosomal replication initiator protein DnaA] [protein_id=WP_000660615.1] [location=197..1558] [gbkey=CDS]"

# headerSplit =  re.split("_", header)
# gene = "_".join(headerSplit[3:5])
# print(gene)

  
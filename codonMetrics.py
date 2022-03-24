#!/usr/bin/env python3
# Name: Aishwarya Basude (abasude)

import sys # allows FastQ to read from STDIN 
import re # regular expressions
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


class FileReader:
    def readFasta(self, file):
        ''' 
        Parses fastA file with header as key, sequence as value. 
        ex: {'YPK_0001': 'atcatc', 'YPK_0002': 'ctgatc', 'YPK_0005': 'atcctgctg'}
        keep in mind every sequence should start with atg and end with stop codon
        '''
        genomeDict = {}

        # parses fasta 
        for line in file:
            line = line.strip() # remove whitespace
            # if header, place as key
            if line[0] == ">":
                # works for seperate formatted file
                # header = re.split(" ", line[1:], 1)
                # gene = header[0]
                # genomeDict[gene] = ""

                # for header in ncbi database cds file
                headerSplit =  re.split("_", line)
                gene = "_".join(headerSplit[3:5])
                genomeDict[gene] = ""

            # if sequence, place as value
            else:
                data = line
                genomeDict[gene] += data

        # special case because of weird output (comment out otherwise)
        # removes invalid genes
        remove = []
        for gene in genomeDict:
            if gene[0] != "W":
                remove.append(gene)

        for gene in remove:
            genomeDict.pop(gene)

        # print(genomeDict)
        return genomeDict
    
    
    def readCpmTpmProkSeq(self, file):
        '''
        Parse cpm tpm data into dictionary with gene name as key and treatment and control tpm as value.
        ex: {'YPK_0001': ['82.1296', '40.1565'], 'YPK_0002': ['141.2988', '197.3905'], 'YPK_0003': ['0.0', '0.0']
        {gene: [sampletpm, controltpm]}
        '''
        inFile = open(file, 'r')
        lines = inFile.readlines()[1:]

        expressionDict = {}
        
        for read in lines:
            splits = re.split(",", read)
            # [sampletpm, controltpm]
            expressionDict[splits[1]] = [splits[13], splits[16]]

        # print(expressionDict)
        return expressionDict
    

    def readCpmTpmSalmon(self, controlTpmfile, treatmentTpmFile):
        '''
        Parse cpm tpm data into dictionary with gene name as key and treatment and control tpm as value.
        ex: {'YPK_0001': ['82.1296', '40.1565'], 'YPK_0002': ['141.2988', '197.3905'], 'YPK_0003': ['0.0', '0.0']
        {gene: [sampletpm, controltpm]}
        '''

        # for treatment file
        treatmentFile = open(treatmentTpmFile, 'r')
        readStressed = treatmentFile.readlines()[1:]

        # for control file
        controlTpmfile = open(controlTpmfile, 'r')
        readControl = controlTpmfile.readlines()[1:]

        expressionDict = {} 

        # add stressed tpm to dict
        for read in readStressed:
            splits = re.split(" +|\t", read)
            header = splits[0]
            headerSplit =  re.split("_", header)
            gene = "_".join(headerSplit[3:5])
            expressionDict[gene] = splits[3]

        # add control tpm to dict checking if key already there
        for read in readControl:
            splits = re.split(" +|\t", read)
            header = splits[0]
            headerSplit =  re.split("_", header)
            gene = "_".join(headerSplit[3:5])
            if gene in expressionDict:
                expressionDict[gene] = [expressionDict[gene], splits[3]]

        # special case because of weird output, removes invalid "genes"
        remove = []
        for gene in expressionDict:
            if gene[0] != "W":
                remove.append(gene)

        for gene in remove:
            expressionDict.pop(gene)

        # print(expressionDict)
        return expressionDict



class CodonMetrics:
    def __init__(self, genomeDict, expressionDict):
        """
        Initialize genome dictionary and expression dictionary from FileReader class.
        """
        self.genomeDict = genomeDict # genes and seqs dictionary
        self.expressionDict = expressionDict # genes and sample and control tpm dictionary
    

    def makeCodonDictDict(self):
        '''
        Count codons for each sequence and store in dictionary.
        Dictionary of dictionaries where first key is gene and each value is a dictionary for codon and normalized codon count.
        Ex: {'YPK_0001': {'ATC': 1.0}, 'YPK_0002': {'CTG': 0.5, 'ATC': 0.5}, 'YPK_0005': {'ATC': 0.3333333333333333, 'CTG': 0.6666666666666666}}
        '''
        codonDictDict = {}
        # iterate through genes
        for gene, sequence in self.genomeDict.items():
            codonDict = {} # hold codon count for each gene
            totalCodonCount = 0  # keeps track of number of codons in each seq

            for i in range(0, len(sequence), 3):
                codon = sequence[i:i+3].upper() # create codon

                # makes sure a triplet
                if len(codon) % 3 == 0: 
                    # if the codon is not already in the dictionary, create key
                    if codon not in codonDict:
                        codonDict[codon] = 0

                    # count codons for each gene and total count
                    codonDict[codon] += 1
                    totalCodonCount += 1

            # normalize codon values by dividing by number of codons in each seq
            for codonKey, numCodons in codonDict.items():
                if totalCodonCount != 0:
                    codonDict[codonKey] = numCodons / totalCodonCount
                
            codonDictDict[gene] = codonDict # add key value pair to dict
        
        return codonDictDict


    def weightGenes(self):
        '''
        Create weighted codon count dictionary of dictionaries.
        Ex: {'YPK_0001': {'ATC': [82.1296, 40.1565]}, 'YPK_0002': {'CTG': [70.6494, 98.69525], 'ATC': [70.6494, 98.69525]}}
        where {gene: {aa: [sample tpm weighted codon count, control tpm weighted codon count]}}
        '''   
        codonDictDict = self.makeCodonDictDict() # create dict of normalized codon values for each gene
        weightedCodonDict = {}

        for gene in codonDictDict:
            # makes sure present in both files/ matches them
            if gene in self.expressionDict:
                codonDict = {} # holds codon dict for each gene 

                for codon, count in codonDictDict[gene].items():
                    # weight the codons by gene expression
                    codonDict[codon] = [float(self.expressionDict[gene][0]) * float(count),
                                        float(self.expressionDict[gene][1]) * float(count)]

                    weightedCodonDict[gene] = codonDict # update dictionary

        return weightedCodonDict


    def makeCodonDict(self):
        '''
        Take in dictionary of dictionarys (genes and codon dicts like weighted codon dict)
        Create dictionary of codon expression for stressed and nonstressed prokaryote across all genes. 
        # example: {'ATC': [223.4284, 237.54700000000003], 'CTG': [141.2988, 197.3905]}
        where {codon: [weighted TPM across all genes for treatment, weighted TPM across all genes control]}
        '''
        weightedCodonDict = self.weightGenes()

        codonDict = {} # holds added weighted codon scores across all genes 
        for gene in weightedCodonDict:
            for codon, value in weightedCodonDict[gene].items():

                # add new weighted codon count without gene to dictionary if not there
                if codon not in codonDict:
                    codonDict[codon] = value
                # update the values for the codons in dict if already there
                else:
                    codonDict[codon] = [codonDict[codon][0] + value[0], codonDict[codon][1]+ value[1]] 

        return codonDict


    def groupedBar(self):
        """Plot bar graph comparing codon expression for stressed and nonstressed prokaryote"""
        codonDict = self.makeCodonDict()

        # codon to amino acid values
        codonToAaDict = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'STOP', 'TAG':'STOP',
            'TGC':'C', 'TGT':'C', 'TGA':'STOP', 'TGG':'W'}

        # amino acid to colors on graph (can change)
        aaColorDict = {
            'M': 'green', 'G': 'blue', 'P': 'orange', 'V': 'red', 'L': 'purple',
            'I': 'brown', 'A': 'gray', 'C': 'olive', 'F': 'pink', 'Y': 'cyan',		
            'W': 'maroon',	'H': 'teal', 'K': 'indigo',	'R': 'lime', 'Q': 'magenta', 
            'N': 'deepskyblue', 'E': 'slateblue', 'D': 'coral', 'S': 'mediumvioletred', 
            'T': 'gold', 'STOP': 'black'}		

        # create data structures to plot on bar
        codons = list(codonDict.keys()) # list of codon names+
        counts = list(codonDict.values()) # list of lists codon values for treatment and control
        treatment = [] # list of treatment weighted codon value
        control = [] # list of control weighted codon values

        # seperate treatment and control values
        for count in counts:
            treatment.append(count[0])
            control.append(count[1])

        # plot grouped bar 
        x = np.arange(len(codons))  # the label locations
        width = .4 # the width of the bars

        plt.bar(x - width/2, control, width, label='control', color = 'cyan' )
        plt.bar(x + width/2, treatment, width, label='treatment', color = 'red')

        # labels
        plt.title('Weighted Codon Counts for Streptococcus pneumoniae for Heat Shock vs. Unstressed Conditions')
        plt.xticks(rotation='vertical')
        plt.xticks(x, codons)
        plt.xlabel("Codons")
        # plt.yscale('log') # did not help
        plt.ylabel("Weighted Codon Counts")
        plt.legend(["Unstressed","Stressed"])

        # color code the codons by aa
        for codonInd in range(len(codons)):
            plt.gca().get_xticklabels()[codonInd].set_color(aaColorDict[codonToAaDict[codons[codonInd]]])

        # legend for amino acid colors
        patchList = []
        aninoAcids = list(aaColorDict.keys())
        for aa in aninoAcids:
            patch = mpatches.Patch(color = aaColorDict[aa], label = aa)
            patchList. append(patch)
        
        plt.legend(handles = patchList, loc='upper left', bbox_to_anchor=(1, 1), prop={'size': 6})
        plt.show()

    
    def topTPMDiff(self):
        '''
        Find genes with biggest disparity of TPM.
        Print table of top 5 genes.
        '''
        # holds list of tupes with of gene and difference between tpm values of sample and control
        differenceTupList = []

        # fill tuple list
        for gene, tpmList in self.expressionDict.items():
            if float(tpmList[1]) - float(tpmList[0]) > 0:
                regulation = "downregulated"
            else:
                regulation = "upregulated"
            differenceTupList.append((gene, abs(float(tpmList[1]) - float(tpmList[0])), regulation))

        # sort list of tuples by difference from greatest to least 
        differenceTupList = sorted(differenceTupList, key=lambda x:x[1], reverse = True)
        
        # prints table for 5 genes with biggest TPM disparity
        print("Gene\t\tTPM Difference\tDirection")
        # change number in range if want more
        for i in range(0,5):
            print("{}\t{}\t{}".format(differenceTupList[i][0], differenceTupList[i][1],  differenceTupList[i][2]))


    def geneCodonDisparity(self):
        '''
        Find mean unweighted normalized codon usage for each gene and sees which have high distance from mean.
        Print table with codon, gene, and disparity.
        '''
        codonDictDict = self.makeCodonDictDict()

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

        # find gene with highest disparity in average codon count for each codon across genes
        # Ex. {'ATC': ['YPK_0001', 0.38888888888888895], 'CTG': ['YPK_0005', 0.2777777777777778]}
        # {codon: [gene with highest difference from avg for codon, difference]}
        # initialize with placeholder values
        geneCodonDisparityDict = {}
        for codon in averageCodonDict:
            geneCodonDisparityDict[codon] = ['gene', 0]

        for gene, codonDict in codonDictDict.items():
            for codon in codonDict: 
                # calculates difference between avg
                difference = abs(averageCodonDict[codon] - codonDict[codon])

                # if difference is bigger than what is already in dict
                if difference > geneCodonDisparityDict[codon][1]:
                    # update dict
                    geneCodonDisparityDict[codon] = [gene, difference]

        # print dict as table
        print("\tCodon\tProtein Product\t   Disparity")
        for codon in geneCodonDisparityDict:
            print("\t{}\t{}\t   {}".format(codon, geneCodonDisparityDict[codon][0], geneCodonDisparityDict[codon][1]))


class CommandLine():
    '''
    Based on David Bernick's code BME 160.
    Handle the command line, usage and help requests.
   
    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.
   
    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output')
        
        # creates arguments for the input quality encodings, and the output quality encodings, making the default output phred 33
        self.parser.add_argument('-ps', '--ProkSeq', action = 'store', nargs='?', const = True, default = False, help = 'uses prok seq output format to parse')
        self.parser.add_argument('-ctp','--cmpTmpFileProkseq', action = 'store', help='input file with cpm and tpm with ProkSeq Output')
        self.parser.add_argument('-tc','--tmpControl', action = 'store', help='input file with TPMs for control prokaryote from salmon output')
        self.parser.add_argument('-tt','--tmpTreatment', action = 'store', help='input file with TPMs for stressed prokaryote from salmon output')

        self.parser.add_argument('-gb', '--groupedBar', action = 'store', nargs='?', const = True, default = False, help = 'prints grouped bar graph')
        self.parser.add_argument('-t', '--topTPMdiff', action = 'store', nargs='?', const = True, default = False, help = 'prints table of genes with greatest TPM difference')
        self.parser.add_argument('-gcd', '--geneCodonDisparity', action = 'store', nargs='?', const = True, default = False, help = 'prints table of gene that causes most disparity for each codon')
        
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


def main(inCL=None):
    '''
    Use STDIN to take in a fastQ file and required arguments specifying the input and output quality encoding.
    Call the CommandLine, the FastQReader, and the Converter classes, and pass them data from the command line.
    Print the header, sequence, description, and translated quality lines. 
    ''' 
    # command line options
    if inCL is None:
        myCommandLine = CommandLine() # no input given to command line class
    else :
        myCommandLine = CommandLine(inCL) # parses arguments if input given 

    # inputs from command lines
    genomeFasta = sys.stdin # takes in genome fna file from STDIN
    
    # call the filereader class
    myReader = FileReader() # calls FastQReader class
    genomeDict = myReader.readFasta(genomeFasta) # calls the readFastQ method from the FastQReader (returns list of dictionary of each quartet)

    # if using prokseq output
    if myCommandLine.args.ProkSeq:
        cmpTmpFileProkseq = myCommandLine.args.cmpTmpFileProkseq # takes in genome expression file from command line 
        expressionDict = myReader.readCpmTpmProkSeq(cmpTmpFileProkseq) # create dictionary for cpm values from file

    # if using salmon output
    else:
        tmpControl = myCommandLine.args.tmpControl
        tmpTreatment = myCommandLine.args.tmpTreatment
        expressionDict = myReader.readCpmTpmSalmon(tmpControl,tmpTreatment)

    # run codon class
    codonData = CodonMetrics(genomeDict, expressionDict)

    # prints based on commandline flags 
    # if -gb or --groupedBar used in commandline as flag
    if myCommandLine.args.groupedBar:
        codonData.groupedBar() # prints grouped bar graph

    # if -t or --topTPMdiff used in commandline as flag
    if myCommandLine.args.topTPMdiff:
        codonData.topTPMDiff() # prints table with genes with highest TPM disparity
        print()

    # if -gdc or --geneCodonDisparity used in commandline as flag
    if myCommandLine.args.geneCodonDisparity:
        codonData.geneCodonDisparity() # prints table with genes that have the highest difference from the mean of codon values for each codon across all genes
        print()


if __name__ == "__main__": # starts functions of main
    main()     

# using prok seq output
# run by inputting: codonCMetrics.py <prokaryoteCDSfastA -ps -ctp cmpTmpfile (optional): -gb -t -gcd 
#"C:\Users\basud\OneDrive\Documents\UCSC\Visual Studio 2019\Project\BME230A\Final\codonCMetrics.py" <"C:\Users\basud\OneDrive\Documents\UCSC\Visual Studio 2019\Project\BME230A\Final\yersinia_pseudotuberculosis_CDS.fasta" --cmpTmpFile "C:\Users\basud\OneDrive\Documents\UCSC\Visual Studio 2019\Project\BME230A\Final\yersinia_pseudotuberculosis_CPM_TMP_sample_count.txt"

# using salmon output
# run by inputting: codonCMetrics.py <prokaryoteCDSfastA -tc tmpControl -tt tmpTreatment  (optional): -gb -t -gcd
# "C:\Users\basud\OneDrive\Documents\UCSC\Visual Studio 2019\Project\BME230A\Final\codonMetrics.py" <"C:\Users\basud\OneDrive\Documents\UCSC\Visual Studio 2019\Project\BME230A\Final\streptococcus_pneumoniae_CDS.fna" -tc "C:\Users\basud\OneDrive\Documents\UCSC\Visual Studio 2019\Project\BME230A\Final\streptococcus_pneumoniae_TPM_control.txt" -tt "C:\Users\basud\OneDrive\Documents\UCSC\Visual Studio 2019\Project\BME230A\Final\streptococcus_pneumoniae_TPM_stressed.txt"
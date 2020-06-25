#!/usr/bin/env python3
# Name: David Lenci (dlenci)
# Group Members: none

class OrfFinder() :
    '''
    The OrfFinder class takes in nucleotide sequences as input and outputs genes based on 
    four parameters: what is a start codon, what is a stop codon, what is the minimum length 
    of the gene, and whether or not to only look for largest gene (This code only looks for larges gene
    and assumes ATG as start codon and TAG, TAA, and TGA as stop codons. Minimum gene length can be 
    chosen).
    Example use:
    commandline: python findORFs.py -lG -s 'ATG' -mG 100 <tass2.fa >tass2ORFdata-ATG-100.txt
    output> Reading frame start..stop length
    
    Psuedocode:
    init: 
        Instantiates OrfFinder class and takes input sequence
        seq = fasta sequence
        
    seqCompliment:
        revSeq = seq backwards
        compSeq = revSeq replacing A:T, T:A, G:C, C:G
        returns compSeq
        
    findOrf:
        seq = orgSeq
        codonLst = seq split every nucleotides
        for codon in codonLst:
            look for stop codon
            if find stop codon:
                create save gene start stop and length in geneComp[]
            if find start codon:
                break
        while codon < len(codonLst):
            if find start codon:
                save start codon location
                while codon in range of start codon position and len(codonLst):
                    if find stop codon:
                        save stop codon location
                        break loop
                    if while loop gets to end of list:
                        save stop codon
                    codon + 1
                codon + stop codon location
            else: codon + 1
            save start position, stop position, and length of gene in list
            save gene to list of genes
        
        adjust positions found of genes so they line up with first reading frame
        
        return gene list
    
    findGene:
        Takes sequence and frames and calls findOrf
        adds frames to each gene in list of genes according to input frame
        returns gene list with frames
    
    makeGeneLst:
        calls seqCompliment and makes sequences for each reading frame (6 total)
        calls findGene for each reading frame with given frame
        combines all gene lists for each reading frame
        
        sorts genes by length decreasing
        
        returns complete list of genes
    '''
    
    def __init__(self, nucSequence):
        
        self.seq = nucSequence
        
        
    def seqCompliment(self, orgSeq):
        '''
        Creates reverse compliment of in put sequence and returns it. Required parameter:
        input sequence.
        parameter: sequence
        
        Example:
        input: AAAGGG
        returns: CCCTTT
        '''
        self.seq = orgSeq
        self.revSeq = self.seq[::-1] #Compliment string
        self.compDict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'} #Dictionary of corresponding nucleotides.
        self.compSeq = ('')
        for i in self.revSeq:
            self.compSeq = self.compSeq + self.compDict[i] 
        return self.compSeq
    
    
    def findOrf(self, seq, orf): 
        '''
        Takes sequence and reading frame and returns all genes giving the start nucleotide,
        end nucleotide, and length of gene.
        Required parameters:sequence and the reading frame.
        
        Example:
        input: ATGAAAAAATAG
        returns: [1, 12, 12]
        
        This code assumes start codon is ATG and stop codons are TAG, TAA, TGA.
        Returns all genes no matter length.
        '''
        self.nucSeq = seq
        self.orf = orf
        self.completeComp = list()
        split = 3
        self.codonLst = [self.nucSeq[i:i+3]for i in range(0, len(self.nucSeq), split)] #Breaks nucleatide into codon list
        self.position = 0
        leng = len(self.codonLst)
        
        self.geneComp = list()
        
        for codon in range(self.position, leng): #Dangling start loop.
            if self.codonLst[codon] == 'ATG': #Stops loop if start found.
                break
            self.position = self.position + 1
            if self.codonLst[codon] == 'TAG' or self.codonLst[codon] == 'TAA' or self.codonLst[codon] == 'TGA':
                self.codonLength = self.position * 3
                self.geneComp+=[1, self.codonLength, self.codonLength] #Creates list of gene start stop and length if one found
                self.completeComp.append(self.geneComp) #Adds gene to list of genes
        
        codon = 0
        while codon < leng: #Loop for finding non-dangling start genes.
            self.position = self.position + 1
            if self.codonLst[codon] == 'ATG':
                loopStart = self.position - 1 
                self.geneStart = self.position
                self.geneStop = 0
                self.geneComp = list()

                while loopStart < leng:
                    if self.codonLst[loopStart] == 'TAG' or self.codonLst[loopStart] == 'TAA' or self.codonLst[loopStart] == 'TGA':
                        self.geneStop = loopStart #Records stop and ends loop when stop found.
                        break
                    else: 
                        loopStart = loopStart + 1
                            
                    if loopStart == leng-1: #Handles dangling stop.
                        self.geneStop = leng
                
                self.nucStart = self.geneStart * 3 - 2 #Transfers codon positions to nucleotide positions
                self.nucEnd = (self.geneStop + 1) * 3
                self.geneLength = self.nucEnd - (self.nucStart-1) #Calculates length
                
                self.geneComp+=[self.nucStart, self.nucEnd, self.geneLength] #Creates gene
                self.completeComp.append(self.geneComp) #Adds gene to gene list
                    
                codon = loopStart + 1
                self.position = loopStart + 1   
            else:
                codon = codon + 1
        
        #Adjusts values for sequences not in first reading frame to first reading frame coordinates.
        if self.orf == '+2':
            for gene in self.completeComp:
                gene[0]=gene[0]+1
                gene[1]=gene[1]+1
        
        if self.orf == '+3':
            for gene in self.completeComp:
                gene[0]=gene[0]+2
                gene[1]=gene[1]+2
                
        if self.orf == '-1' or self.orf == '-2' or self.orf == '-3':
            leng = len(self.nucSeq)
            for gene in self.completeComp:
                gene[0] = leng - gene[1] + 1
                gene[1] = gene[0] + gene[2] - 1 
            
        return self.completeComp
    
    
    def findGene(self, seq, orf):
        '''
        Calls findOrf and adds frame to gene list. Returns complete gene list
        with reading frames for each gene.
        Parameters: sequence nd reading frame for sequence.
        
        Example:
        input: ATGAAAAAATAG, +1
        returns: [1, 12, 12, +1]
        '''
        self.nucSeq = seq
        self.orf = orf
        self.geneLst = self.findOrf(self.nucSeq, self.orf) #Creates gene list
        for gene in self.geneLst: #Adds frame to each gene in gene list.
            gene.append(self.orf)
            
        return self.geneLst
    
    
    def makeGeneLst(self):
        '''
        Establish the reading frames for sequence and its compliment. Calls complimentSeq
        and establishes all reading frames. Calls findGene function and sorts list of genes.
        Takes sequence from init function for calls.
        
        returns: [[start, stop, length, frame]]
        '''
        #Establish positive reading frames
        self.orfOne = self.seq
        self.orfTwo = self.seq[1:]
        self.orfThree = self.seq[2:]
        
        self.compSeq = self.seqCompliment(self.seq) #Creates reverse compliment of seq
        
        #Establish reading frames for reverse compliment
        self.cOrfOne = self.compSeq
        self.cOrfTwo = self.compSeq[1:]
        self.cOrfThree = self.compSeq[2:]
        
            
        #All frames
        plusOne = "+1"
        plusTwo = "+2"
        plusThree = "+3"
        minusOne = "-1"
        minusTwo = "-2"
        minusThree = "-3"
        
        #Finds genes in each frame    
        self.genesInOne = self.findGene(self.orfOne, plusOne)
        self.genesInTwo = self.findGene(self.orfTwo, plusTwo)
        self.genesInThree = self.findGene(self.orfThree, plusThree)
            
        self.genesIncOne = self.findGene(self.cOrfOne, minusOne)
        self.genesIncTwo = self.findGene(self.cOrfTwo, minusTwo)
        self.genesIncThree = self.findGene(self.cOrfThree, minusThree)   
        '''Combine gene lists and sort based off of length of genes.'''
         
        self.completeGeneLst = []
        
        #Combines all found genes
        self.completeGeneLst.extend(self.genesInOne)
        self.completeGeneLst.extend(self.genesInTwo)
        self.completeGeneLst.extend(self.genesInThree)
        self.completeGeneLst.extend(self.genesIncOne)
        self.completeGeneLst.extend(self.genesIncTwo) 
        self.completeGeneLst.extend(self.genesIncThree)

        
        #Sorts genes in decreasing order   
        self.completeGeneLst = sorted(self.completeGeneLst, key=lambda x: (x[2], x[0]), reverse=True)
        
        
        return self.completeGeneLst
    
class CommandLine():
   
    """
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.
    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    """
    def __init__(self, inOpts=None):
        """
         CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        :param inOpts: 
        """
        import argparse
        self.parser = argparse.ArgumentParser(
            description='Finds the largest ORF in a DNA sequence',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store_true', help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= range(0, 2000),
            action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action='append', nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)
            


import sequenceAnalysis

def main(myCommandLine=None):
    '''
    Reads in a fasta file and outputs the ORFs frame, start, stop, and length in specified ourput file.
    If no arguments are give assumes default values: min length of 100, start is ATG, stop is TAG, TAA, TGA
    and longest gene.
    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()
        print (myCommandLine.args)
        
        fastaFile = sequenceAnalysis.FastAreader()
            
        for header, sequence in fastaFile.readFasta():
            print(header)
            orfData = OrfFinder(sequence)
            orfLst = orfData.makeGeneLst()
            filteredList = filter(lambda orf: orf[2] > myCommandLine.args.minGene, orfLst)# Removes orfs lower than min length.
                
            for orfs in filteredList:  
                frame = orfs[3]
                start = orfs[0]
                stop = orfs[1]
                length = orfs[2]
                print('{:s} {:>5d}..{:>5d} {:>5d}'.format(frame, start, stop, length))#Prints out data in specific order

  
 
if __name__ == "__main__":
    main()

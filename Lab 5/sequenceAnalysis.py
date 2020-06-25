#!/usr/bin/env python3
# Name: David Lenci (Cdlenci)
# Group Members: none

import collections 

class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        '''
        Initializes ProteinParam object, and creates a dictionary with keys for every amino acid 
        and values for the number of that acid key in the input sequence.
        Paramater:Protein sequence string.
        '''
        upperProtein = protein.upper() # Makes sure protein string is upper case.
        self.newdict = {}
        for aa in self.aa2mw.keys(): # Creates new dictionary using keys from aa2mW
            self.newdict[aa] = 0 # Assigns keys in new dictionary to 0
        for x in upperProtein: # Iterates through protein string.
            if x in self.newdict:
                self.newdict[x] = self.newdict[x] + 1 # Increases value in dictionary everytime 'key' occurs in string.
            else:pass

    def aaCount (self):
        '''
        Sums all amino acid values to get total aa count.
        Returns that sum.
        '''
        self.aaSum = sum(self.newdict.values()) # Sums all values in aa comp dictionary.
        
        return self.aaSum
        
    def pI (self):
        '''
        Iterates over all values of pH from 0-14 to the hundreths place.
        Determines the pH in that range that is closest to net neutral charge.
        Uses charge function to determine the charge of protein sequence at each pH value.
        '''
        self.chargeOld = 10000000
        self.pHBest = 0
        for i in range(1401): # By making range 0-1400 and then dividing by 100 we get pH values to the hundreths place.
            pH = i/100
            self.chargeNew = self.charge(pH)
            if abs(self.chargeNew) < abs(self.chargeOld): # Finds charge closest to 0.
                self.chargeOld = self.chargeNew
                self.pHBest = pH
        return self.pHBest

    def aaComposition (self) :
        '''
        Returns dinctionary with amino acid composition of protein sequence.
        '''
        self.aaComp = {}
        for aa in self.newdict.keys():
            self.aaComp[aa] = self.newdict[aa]
        else:pass
        return self.aaComp # Returns dictionary made in init.
    

    def charge (self, pHvalue):
        '''
        Returns the charge of protein sequence at a given pH value.
        Parameters:pH input.
        '''
        self.pH = pHvalue
        self.totalNeg = 0
        self.totalPos = 0
        for aa in self.aa2chargePos.keys(): # Iterates through all AAs with pos charge.
            self.chargePos = (self.newdict[aa]*((10**self.aa2chargePos[aa])/((10**self.aa2chargePos[aa])+10**self.pH))) 
            self.totalPos = self.chargePos + self.totalPos # Returns total pos charge by relevant AAs.
        
        for aa in self.aa2chargeNeg.keys(): # Iterates thorugh all AAs with neg charge.
            self.chargeNeg = self.newdict[aa]*((10**self.pH)/((10**self.aa2chargeNeg[aa])+10**self.pH))
            self.totalNeg = self.chargeNeg + self.totalNeg # Returns total neg charge by relevant AAs.
                                               
        self.totalPos = (10**self.aaNterm/(10**self.aaNterm+10**self.pH)) + self.totalPos # Adds N-terminus to pos charge.
        self.totalNeg = (10**self.pH/(10**self.aaCterm+10**self.pH)) + self.totalNeg # Adds C-terminus to neg charge.
        self.netCharge = self.totalPos - self.totalNeg # Calculates and returns net charge at input pH.
        return self.netCharge

    def molarExtinction (self):
        '''
        Return molar extinction coefficient for given protein sequnce.
        '''
        self.mExtinc = self.newdict['Y']*self.aa2abs280['Y']+self.newdict['W']*self.aa2abs280['W']+self.newdict['C']*self.aa2abs280['C']
        return self.mExtinc # Returns extinction coefficient using given values for Y, W, and C.

    def massExtinction (self):
        '''
        Converts molar extinction coefficient to mass extinction coefficient.
        '''
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0 # Returns conversion of extinction coefficient.

    def molecularWeight (self):
        '''
        Calculates and returns molecular weight of protein sequence.
        '''
        self.mWeight = {}
        for aa in self.aa2mw.keys():
            self.mWeight[aa] = self.newdict[aa]*self.aa2mw[aa] # Calculates molecular weight using given dictionary.
            
        self.mwSum = sum(self.mWeight.values())
        waterCount = self.aaSum - 1
        waterWeight = waterCount*self.mwH2O # Subtracts weight water for every peptide bond.
        self.finalWeight = self.mwSum - waterWeight
        return self.finalWeight

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    


class NucParams:
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        self.nucSeq = inString
        self.upperNucSeq = self.nucSeq.upper()
        self.newNucSeq = self.upperNucSeq.replace('T','U')
        self.pureNucSeq = self.newNucSeq.replace(' ','')
        self.rnaNucSeq = self.pureNucSeq
        
        self.codeCount = {}
        split = 3
        self.codonLst = [self.rnaNucSeq[i:i+3]for i in range(0, len(self.rnaNucSeq), split)]
        for code in self.rnaCodonTable.keys():
            self.codeCount[code] = 0
        
        self.aaComp = {}
        allAA = 'AGMSCHNTDIPVEKQWFLRY-'
        sortedAA = ''.join(sorted(allAA))
        for aa in allAA:
            self.aaComp[aa]=0
            
         
    def addSequence (self, inSeq):
        pass
    
    def aaComposition(self):
        for aa in self.aaComp.keys():
            for code in self.rnaCodonTable.keys():
                if aa == self.rnaCodonTable[code]:
                    self.aaComp[aa]=self.aaComp[aa]+self.codeCount[code]
                else:pass
        
        return self.aaComp
    
    def nucComposition(self):
        self.nucComp = {}
        self.NucSeq = self.nucSeq.upper()
        validNuc = 'ACGTUN'
        
        for n in validNuc:
            self.nucComp[n]=0
        for n in self.NucSeq:
            self.nucComp[n]=self.nucComp[n]+1
            
        return self.nucComp
    
    def codonComposition(self):
        for code in self.codonLst:
            if code in self.codeCount:
                self.codeCount[code] = self.codeCount[code] + 1
            else: pass
            
        return self.codeCount
    
    def nucCount(self):
        self.validSum = sum(self.nucComp.values())-self.nucComp['N']
        return self.validSum
    
    

import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
        

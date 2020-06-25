class tRNA():
    '''
    Takes a list of tRNA sequences and finds all essential kmers of each sequence based 
    off of the other tRNA sequences input. Does this by first finding the power set of each sequence 
    then reducing this to only unique sequences. This set is then reduced to only essential 
    sequences.
    
    Example:
    input:fasta file of tRNA:
          > tRNA | Glu | ∃UC | Bos taurus | mitochondrial
          -GUUCUUGU"GUUGAA---UG---ACAACLAPGGUUU∃UCAUAPCAUUA-------------------G-U?AUGGUPAG"UUCCAUGUAAGAAUACCA 
         ...
    return:
          tRNA|Glu|∃UC|Bostaurus|mitochondrial
          GUUCUUGU"GUUGAAUGACAACLAPGGUUU∃UCAUAPCAUUAGU?AUGGUPAG"UUCCAUGUAAGAAUACCA
          GUUCUU
          .UUCUUG
          ...CUUGU
          ......GU"GUU
          ........"GUUG
          ..........UUGAA
          ............GAAUG
          ..............AUGAC
          ................GACAAC
          ....................ACL
          ......................LA
          .......................APGG
          ........................PGGU
          ...........................UUU∃
          ..............................∃UC
          ................................CAUAP
          .................................AUAPC
          ..................................UAPCA
          ...................................APCAU
          .......................................UUAGU
          ..........................................GU?
          ...........................................U?A
          ............................................?AU
          ..............................................UGGUP
          ................................................GUPA
          ...................................................AG"
          ....................................................."UUCCA
          ......................................................UUCCAU
          .........................................................CAUG
          ............................................................GUAAGA
          .............................................................UAAGAA 
    
    ...
    
    tRNA are output in alphabetical order mimicking the above output and dashes 
    and underscores are removed from output sequence. 
    '''
    def __init__(self, tRNAs):
        '''
        Instantiates tRNA class, and calls the power sets function to create a list
        of tRNA objects where each object contains a header and sequence and associated
        power set. 
        Requires an input that is a list of tRNAs to analyzed, which is input when
        class is instantiated. 
        '''
        self.tRNAlist = tRNAs
        
        for tRNA in self.tRNAlist:
            tRNA[1]=tRNA[1].replace('-','') #remove '-' and '_' from sequences
            tRNA[1]=tRNA[1].replace('_','')
        
        self.powerSets = []
        
        for tRNA in self.tRNAlist: #Create powerset for each tRNA object
            self.pSet = self.powerSet(tRNA[1])
            self.powerSets.append(self.pSet)
        
    def powerSet(self, seq):
        '''
        Given a list of sequences finds the powerset for each sequence.
        Example:
        input: [AGTC, CTGA]
        returns: [{AGT, AG, A, GTC, GT, G, TC, T, C, ATGC}, {CTG, CT, C, TGA, TG, GA, G, A, CTGA}]
        '''
        self.powSet = set()
        length = len(seq)
        for position in range(length):
            for i in range(length):
                if (position + i) < (length):
                    self.powSet.add(seq[position:position + i + 1])        
        return self.powSet
    
    def uniqify(self, alist):
        '''
        Finds the unique items in each powerset in a list of  powersets.
        
        Takes a list of powersets and finds all unique items in each powerset by 
        finding the union of all other powersets for each powerset and subtracting 
        that powerset from the intersection of that powerset with the union of all 
        other powersets.
        
        input: [{AGT, AG, A, GTC, GT, G, TC, T, C, ATGC}, {CTG, CT, C, TGA, TG, GA, G, A, CTGA}]
        returns: [{AGT, AG, GTC, GT, TC, AGTC}, {CTG, CT, TGA, TG, GA, CTGA}]
        '''
        unique = alist
        self.uniqueList = []
        for aSet in unique:
            fullUnion = set()
            uniqueCopy = unique.copy() #Creates copy of powerset list
            uniqueCopy.remove(aSet) #Removes set being uniqified
            
            for pSet in uniqueCopy:
                fullUnion = fullUnion.union(pSet) #Creates union of all other power sets
            intersect = aSet.intersection(fullUnion) #Finds intersection between union and each set.
            uniques = aSet.difference(intersect) #Finds difference between each set and intersection.
            
            self.uniqueList.append(uniques)
            
        return self.uniqueList
        
    def makeEssential(self):
        '''
        Finds all nonessential items in a unique set and subtracts those 
        items from said set creating a set of essential items. Calls on 
        uniqify to create list of unique sets to find essentials from. 
        Essential items are created by extending the length of the unique 
        item in both directions and removing all these items from the unique
        item list.
        
        Example:
        Input: [{AGT, AG, GTC, GT, TC, AGTC}, {CTG, CT, TGA, TG, GA, CTGA}]
        returns: [{AG, GT, TC}, {CT, TG, GA}]
        '''
        
        uniques = self.uniqify(self.powerSets) #Creates list of unique sets
        seqID = 0 #Creates ID for original sequence
        self.essList =[]
        for uniqueSet in uniques:
            seq = self.tRNAlist[seqID][1] #Finds original sequence
            nonEss = set()
            for kmer in uniqueSet: #For each kmer extends kmer to right and left
                leng = len(kmer)
                position = seq.find(kmer)
                end = position + leng
                extend1 = seq[position-1:end]
                extend2 = seq[position:end+1]
                nonEss.add(extend1)
                nonEss.add(extend2)
                
            essentials = uniqueSet.difference(nonEss) #Subtracts nonessential list from unique list to give essential list
            self.essList.append(essentials)
            seqID+=1
        return self.essList

def main():
    import sequenceAnalysis as sa
    
    myReader = sa.FastAreader()
    
    tRNAlst = []
    for head, seq in myReader.readFasta() :
        tRNAOb = [head.lstrip(), seq]
        tRNAlst.append(tRNAOb)
        
    tRNAlst.sort(key=lambda x:x[0]) #tRNA list sorted by tRNA name
    
    mytRNA = tRNA(tRNAlst)
    finalList = mytRNA.makeEssential() #Creates essential tRNA list
    for item in tRNAlst:
        print(item[0])
        print(item[1])
        seq = item[1]
        for fset in finalList: #Prints essential items in readable output
            kmerlist = []
            for kmer in fset:
                pos = seq.find(kmer)
                alist = [kmer, pos]
                kmerlist.append(alist)
            kmerlist.sort(key=lambda x:x[1])
            for item in kmerlist:
                if item[1] != -1:
                    output = '.'*item[1] + item[0]
                    print(output)

if __name__ == "__main__":
    main() 

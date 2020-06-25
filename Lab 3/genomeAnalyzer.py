def main ():
    import sequenceAnalysis as sa
    
    myReader = sa.FastAreader() # make sure to change this to use stdin C:\\Users\\dlenc\\Downloads\\testGenome.fa
    totalSeq=''
    for head, seq in myReader.readFasta() :
        totalSeq=totalSeq+seq
    
    seq=totalSeq
    myNuc = sa.NucParams(seq)
    codonComp = myNuc.codonComposition()
    aaComp = myNuc.aaComposition()
    rnaCodonTable = {
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    
    nucs = sorted(rnaCodonTable.items(), key=lambda value: (value[1], value[0]))
    
    nucComp = myNuc.nucComposition()
    nucSum = myNuc.nucCount()
    gTot=nucComp['G']
    cTot=nucComp['C']
    gcComp = gTot+cTot
    gcPer = (gcComp/nucSum)*100
    nucMb = nucSum/1000000
    
    print("Sequence length = {:.2f} Mb".format(nucMb))
    print('')
    print("GC content = {:.1f}%".format(gcPer))
    print('')
        
    # calculate relative codon usage for each codon and print
    for nucI in nucs:
        nuc = nucI[0]
        aa = nucI[1]
        val = codonComp[nucI[0]]/aaComp[nucI[1]]
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(nuc, aa, val*100, codonComp[nucI[0]]))

if __name__ == "__main__":
    main()

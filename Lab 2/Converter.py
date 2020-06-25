#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python3 
# Name: David Lenci (dlenci) 
# Group Members: none

'''
Program docstring goes here
'''
short_AA = {
            'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
            }


long_AA = {value:key for key,value in short_AA.items()}

RNA_codon_table = {
# Second Base
# U             C             A             G
#U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
#C 
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
#A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
#G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}


dnaCodonTable = {key.replace('U','T'):value for key, value in RNA_codon_table.items()}


def main():
    '''
    Takes either single code for AA, triple code for AA, or codon sequence in either 
    DNA or RNA nucleotides as input, and converts accordingly:
    Single code AA is converted to a triple code AA
    triple code AA is converted to a single code AA
    DNA or RNA codon is converted to a triple code AA
    '''
    code = input("Enter a sequence amino acid value:")
    upperCode = code.upper()
    lng=int(upperCode.count(''))
    codeParseB='same'
    '''
    'if else' parsing statement that passes input through each dictionary to determine 
    if there are any corresponding key values. If there isn't any returns 'Unknown'.
    '''
    if lng == 2:
        print(upperCode+' = '+long_AA.get(upperCode,'Unknown'))
        #Assumes that if input was one letter then input is single code AA.
    if lng == 4:
        codeParseA = RNA_codon_table.get(upperCode,'Unknown')
        #Since there are three possible dictionaries corresponding to a 3 letter 
        #input, code passes that input through each dictionary before printing 'Unknown.'
        if codeParseA == 'Unknown':
            codeParseB=dnaCodonTable.get(upperCode,'Unknown')
        else:
            print(upperCode+' = '+codeParseA)
            
        if codeParseB == 'Unknown':
            codeParseC=short_AA.get(upperCode,'Unknown')
            print(upperCode+' = '+codeParseC)
        else:
            if codeParseB == 'same': #If codeParseB was never assigned then code doesn't need to check its value.
                pass
            else:
                print(upperCode+' = '+codeParseB)
    
main()

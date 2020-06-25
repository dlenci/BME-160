#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python3 
# Name: David Lenci(dlenci) 
# Group Members: None

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example: 
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''

class DNAstring (str):

    def purify(self):
        ''' Return an upcased version of the string, collapsing a single run of Ns. Assuming there is only one chunk of Ns'''
        nUpper = self.upper()
        nPosition = nUpper.find('N')
        nLength = str(nUpper.count('N')) 
        nReplace = nUpper.replace(nUpper[nPosition],'{'+nLength+'}', 1) #Remove the first N in chain of Ns, and replace with N count.
        nLess = nReplace.replace('N','') #Remove remaining Ns
        return (nLess)
    
def main():
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?')
    thisDNA = DNAstring (data)
    pureData = thisDNA.purify()
    print (pureData)
    
main()

#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python3 
# Name: David Lenci (dlenci) 
# Group Members: none

'''
Takes FASTQ sequence name and breaks it down into 
its corresponding components and identifies those components.
'''

class FastqString (str):
    ''' 
    Splits and identifies components of fastq name 
    with one string parameter required'''
    def parse(self):
        '''
        Splits string wherever there is a ':', and removes '@'. 
        Each aspect of the newly created list is then identified, and printed to output.
        '''
        removeAt = self.replace('@','')
        cPosition = removeAt.split(':')
        print('Instrument = '+cPosition[0])
        print('Run ID = '+cPosition[1])
        print('Flow cell ID = '+cPosition[2])
        print('Flow cell lane = '+cPosition[3])
        print('Tile Number = '+cPosition[4])
        print('X-coord = '+cPosition[5])
        print('Y-coord = '+cPosition[6])

def main():
    ''' Function docstring goes here'''
    fastQ = input('FASTQ file?')
    thisFastQ = FastqString (fastQ)
    parseFastQ = thisFastQ.parse()
    
   
main()

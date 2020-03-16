# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 12:32:10 2020

@author: Dani
"""

# packages
import os
from optparse import OptionParser
 


# Release information
__version__ = '0.2'
_scriptname = 'biospa_header'
_verdata = 'Feb 2020'
_devflag = True



# Option parser

parser = OptionParser()
parser.add_option("-i", "--input", 
					dest = "inputfile", 
					metavar = "INPUT",
                  	help = "input file to copy pattern")
parser.add_option("-o", "--output", 
					dest = "outfolder",  
					metavar = "OUTPUT", 
                  	help = "output folder where files will be saved",
                    default = 'Output')

(options, args) = parser.parse_args()




# Program Header
print('\n====================================================\n')
print(_scriptname + ' script, v' + __version__ , _verdata + 
	'\n =-= by Daniel Martinez =-=')
if(_devflag):
    print('\nWARNING! THIS IS JUST A DEVELOPMENT SUBRELEASE.' +
          '\nUSE IT AT YOUR OWN RISK!\n')
print('\n====================================================\n')


# create folder
os.mkdir(str(options.outfolder))

path = os.getcwd()

outdir = path + '\\' + str(options.outfolder)


# open pattern file and store times as variable
source = open(str(options.inputfile), 'rb')
f1 = source.readlines()
times = f1[2]
source.close()

# list all the files with txt extension
files = []
for file in os.listdir(path):
    if file.endswith(".txt"):
        files.append(file)

for file in files:
    dest = open(file, 'rb')
    f2 = dest.readlines()
    new_file = outdir + '\\' + str(file)
    f2[2] = times
    # print(new_file)
    dest2 = open(new_file, 'wb')
    for line in f2:
        dest2.write(line)
    dest.close()
    dest2.close()


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
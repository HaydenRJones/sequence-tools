# flyeLogParser.py
# HJ - 21/10/24
# Parse log output from flye and colate into a single excel file
# Spit out contigs with a single assembly into fasta files too!

# v0.2 - 15/02/25
# - added argument handling
# - general tidy up

# v0.3 - 08/03/25
# - added logging for loci stats (in prep for MAR.py)

#  todo:
# - plotting
# - better file handling etc.

# Imports

import argparse
#import os
#import subprocess

import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns

from Bio import SeqIO

# Functions

def getRef(refFile):
    
    records = list(SeqIO.parse(refFile, 'fasta'))
    
    refList = []
    
    for record in records:
        refList.append(str(record.id))
    
    return(refList)

def openFile(sample, gene):
    
    try:
    
        file = pd.read_csv(f'./flye_assembly/{sample}/assembly/{gene}/flye/assembly_info.txt', sep = '\t', header = 0)
        nContigs = len(file)
    
    except:
        
        nContigs = 0
    
    return(nContigs)

# Main

if __name__ == '__main__': 
    
    parser = argparse.ArgumentParser(
        prog = 'flyeLogParser.py',
        description = 'Parse log output from flye and colate into a single excel file\nAlso get loci with a single contig across all samples and write these to a fasta file')        
    parser.add_argument('-r', '--ref', 
                        help = 'Fasta ref file')
    parser.add_argument('samples', nargs = '+', 
                        help = 'List of samples to parse')
    args = parser.parse_args()
    
    # Loop

    refList = getRef(args.ref)

    logFile = np.zeros((len(refList), len(args.samples)))
    logFile = pd.DataFrame(logFile, columns = args.samples, index = refList)

    for sample in args.samples:
        
        for ref in refList:
            
            logFile.loc[ref, sample] = openFile(sample, ref)

    oneContigLoci = logFile[(logFile == 1).all(axis  = 1)]
    oneContigLoci = oneContigLoci.index.values.tolist()
    
    # Log files
    logFile.to_csv('./flye_assembly/flye_contig_count.tsv', sep = '\t')
    pd.DataFrame(oneContigLoci).to_csv('./flye_assembly/flye_one_contig_loci.tsv', sep = '\t', header = False, index = False)
    logFile[~(logFile.eq(1).all(axis=1))].to_csv('./flye_assembly/flye_multi_contig_loci.tsv', sep = '\t')
    
    for sample in args.samples:
      
        temp = []
      
        for loci in oneContigLoci:
      
            sequence = list(SeqIO.parse(f'./flye_assembly/{sample}/assembly/{loci}/flye/assembly.fasta', 'fasta'))
            temp.append([f'{loci}', f'{sequence[0].seq}'])

        #SeqIO.write(temp, f'{sample}_contigs.fasta', 'fasta')
        file = open(f'./flye_assembly/{sample}_contigs.fasta', 'w')
        for i in range(len(temp)):
            file.write(f'>{sample}_{temp[i][0]}\n{temp[i][1]}\n')

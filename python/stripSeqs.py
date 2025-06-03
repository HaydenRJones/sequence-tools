# stripSeqs.py
# parse meshclust3 output into induvidual fasta files
# HJ - 31/10/23

import pandas as pd
from Bio import SeqIO

def stripSeqs(barcode, clusterDist)

        colnames = ['cluster', 'index', 'identity', 'score']
        clusters = pd.read_csv(f'./clusters/{clusterDist}_{barcodes}.txt', sep = '\t', header = None, names = colnames)
        clusters = clusters[clusters.groupby('cluster')['cluster'].transform('size') > 1]
        
        with open(f'./aligned/{barcode}.fas') as fasta_file:
            identifiers = []
            sequences = []
            for record in SeqIO.parse(fasta_file, 'fasta'):
                identifiers.append(record.id)
                sequences.append(str(record.seq))
     
        sequences = {'index': identifiers, 'sequence': sequences}
        sequences  = pd.DataFrame(sequences)
        
        for x in range(1, clusters['cluster'].nunique() + 1):
                
            clusterSeqs = pd.DataFrame(columns = ['index', 'sequence'])
        
            workingData = clusters.loc[clusters['cluster'] == x]
        
            for j in range(0, len(workingData)):
                target = workingData['index'].iloc[j].strip('>')       
                seq = pd.DataFrame(sequences.loc[sequences['index'] == target])
        
                if not seq.empty:
                    
                    clusterSeqs = pd.concat([clusterSeqs, seq], ignore_index = True , keys=['index', 'sequence'])
        
                else:
                    print(f'Index \'{target}\' not found in sequences DataFrame')
                    
                
                f = open(f'./clusters/{barcode}_{x}.fasta', 'w')
                
                for k in range(0, len(clusterSeqs)):
                    
                    name = clusterSeqs['index'].iloc[k]
                    seq = clusterSeqs['sequence'].iloc[k]
                    
                    f.write(f'>{name}\n{seq}\n')
                    
                f.close()
                
            print(f'Clusters stripped for {barcode}, cluster {x}')
    
    return()        

#########

barcodes = []
clusterDist = 99

for i in range(0,len(barcodes))
    
    stripSeqs(barcodes[i], clusterDist)

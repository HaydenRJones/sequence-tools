# nQuireDet.py
# HJ - 07/04/26
# v0.3 
# Updated from nQuireDet.py (07/05/25)
# - Removed code that splits the input reference file based on ploidy
# - Removed requirement for a reference file
# Updated from 0.2
# - Check for and flag loci without variation. This shows up as a freemodel with a score of zero or -nan

import pandas as pd
import sys
#import argparse

if __name__ == '__main__': 
    
    sample = sys.argv[1]
    
    df = pd.read_csv(f'{sample}_nQuire_stats.tsv', sep = '\t')

    ploidy = ['dip', 'tri', 'tet']

    # name, det, dLogLK
    det = []
    
    for i in range(0, len(df)):
    
        values = list((df.iloc[i])[5:8])
        idx = values.index(min(values))
        
        #this sucks!    
        name = (df.iloc[i][0]).split('/')[-1]#((df.iloc[i][0]).split('/')[-1]).split('_')[:-1]
        #name = '_'.join(name)
        
        freemodel = (df.iloc[i][2])
        
        if freemodel == 0 or freemodel == '-nan' or freemodel == 'nan':
            det.append([name,'NA','NA'])
            continue
        
        det.append([name,
                    ploidy[idx],
                    float(min(values))])

    det_df = pd.DataFrame(det, columns = ['name', 'det', 'dLogLK'])
    det_df.to_csv(f'{sample}_det.tsv', sep = '\t', index = False)
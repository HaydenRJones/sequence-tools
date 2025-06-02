# resolveMultiContigs.py
# HJ - 06/04/25
# Take in the multi contig list from flyeLogParser.py and attempt to resolve loci with multiple contigs

# Imports

import pandas as pd
import numpy as np
import mappy as mp
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolours

import argparse

from Bio import SeqIO

# Funcs

def make_Alignment(directory, sample, loci, ref_list):
    
    index = mp.Aligner(f'{directory}/{sample}/assembly/{loci}/flye/assembly.fasta', preset='splice')
    hit_temp = []
    cigar_array = []
    
    qer_length = 0
    
    for name, seq in ref_list:
    #for name, seq, qual in mp.fastx_read(f'{reference}'):
        for hit in index.map(seq):
            
            #print(hit.cigar_str)
            
            if name == loci and qer_length == 0: qer_length = len(seq)
            
            # Manually add in clipping cigar strings infered using the start and end of the query mapping.
            # It looks like there might be a bug in mappy as this should be reported but isn't. 
            # If the bug is fixed and clipping is reported properly !this will break!
            cg = hit.cigar
            cg.insert(0, [hit.q_st, 4])
            #cg.append([hit.q_en, 4])    # !This is wrong! Something weird is going on need
            cg.append([qer_length - hit.q_en, 4])
            
            cigar_array.append(hit.cigar)
            
                        
            hit_temp.append([hit.ctg, 
                             name, 
                             hit.r_st, 
                             hit.r_en, 
                             hit.r_en - hit.r_st, 
                             hit.q_st, 
                             hit.q_en, 
                             hit.q_en - hit.q_st, 
                             hit.NM, 
                             (hit.q_en - hit.q_st) / len(seq), 
                             hit.mapq,
                             cg])
            
    hit_table = pd.DataFrame(hit_temp, columns = ['CONTIG', 
                                                  'QER_NAME', 
                                                  'REF_START', 
                                                  'REF_END' , 
                                                  'MAP_LEN', 
                                                  'QER_START', 
                                                  'QER_END', 
                                                  'QER_MAP', 
                                                  'QER_MISMATCH',  
                                                  'LEN_MAP', 
                                                  'MAPQ',
                                                  'CIGAR'])
    
    #plot_Cigar(cigar_array, list(hit_table['CONTIG']), qer_length, loci)
    
    return(hit_table, qer_length)

def find_Containation(hit_table, contigs, loci_name):
    
    bad_contigs = []
    
    #contigs = list(assemblyInfo['#seq_name'])
    
    hits = list(hit_table['CONTIG'])
    
    unmapped_hits = set(contigs) - set(hits)
    for x in unmapped_hits:
        bad_contigs.append([x, 'unmapped'])
    
    for x in hit_table['QER_NAME']:
        if x != loci_name:
            bad_contigs.append([x, 'wrong_loci'])
    
    # Take our hit table and find if any contigs had nothing map to it -> discard these. They might be something but if paftol doesn't map nothing we can do about that
    # Also check the names of the queries and if these match the loci we're checking -> if a different paftol maps than expected that suggestes contamination too
    
    # Is also going to need a list of possible contigs to check -> get this from our assemblyinfo.txt, mapped contigs from the hit table
    # Should we also consider the read depth reported during assembly? we might assume this is low(er) than true assemblies? 
    
    return(bad_contigs)

# For these next two maybe we could start by computing pairwise overlap. -> from this we could try clustering it, but might be super overkill
# Paralogs would have signficant overlap 
# Disjunct assembilies wont

def find_Partial(hit_table, thresh):   
    
    partial_contigs = []
    
    for i in range(0, len(hit_table)):
        if hit_table.iloc[i]['LEN_MAP'] <= thresh:
            partial_contigs.append(hit_table.iloc[i]['CONTIG'])            
    
    return(partial_contigs)

def pick_Paralog(hit_table):

        

    return()

def plot_Cigar(cigar_arr, contigs, query_length, loci):
    
    # Make an empty array that is the length of our paftol reference
    hit_map = np.zeros((len(cigar_arr), query_length), dtype = int)
    
    for i in range(len(hit_map)):
        
        start_pos = 0
        
        for string in cigar_arr[i]:
            
            # If we are dealing with a clip or match add that to our hit map. In total these should add to the length of the paftol reference so simply add as is.
            if string[1] == 0:
                hit_map[i][start_pos:start_pos + string[0]] = 1
                start_pos += string[0]
            if string[1] == 4:
                  hit_map[i][start_pos:start_pos +string[0]] = 0
                  start_pos += string[0]  
            
            # In the case of a splice (infered intron) we don't really have room for that in our map, so what we do is replace the last position in the hit map
            # This means or final hit map isn't strictly true as all exons are shown as 1bp shorter, but it lets us visualse things equally.
            if string[1] == 3:
                hit_map[i][start_pos] = 2
                start_pos += 1
    
            # CIGAR flags 1 ans 2 are indels in the reference or query. These are just skipped for the sake of simplicty
            # We could probably look at adding this if it's really important but for now it's too much work
            # If implemented it would be basically a copy paste of the splice check.
    
    # Red Green Purple
    cmap = mcolours.ListedColormap(['#914343', '#28a13f', '#6d47ad'])
                
    plt.figure(dpi=300)
    sns.heatmap(hit_map, 
                vmin = 0, 
                vmax = 2, 
                cmap = cmap, 
                cbar = False,
                yticklabels = contigs).set_title(f'{sample}/{loci}')
    # Add a horizontal line between each hit to break things up a bit
    for x in  range(1, len(hit_map)):
        plt.axhline(x, c = 'black', lw = 0.5)    
    plt.show()
    
    return()

# Main

if __name__ == '__main__': 
    
    # parser = argparse.ArgumentParser(
    #     prog = 'resolveMultiAssembly.py',
    #     description = '')
        
    # parser.add_argument('-r', '--ref', 
    #                     help = 'Fasta ref file')
    # parser.add_argument('-f', '--flye-multi-contigs', 
    #                     help = 'multi contig file from flyeLogParser')
    # parser.add_argument('-d', '--directory', 
    #                     help = 'assembly directory, defaults to ./flye_assembly/SAMPLES/assembly/LOCI')
    # args = parser.parse_args()
    
    reference = './references/paftolref_Libertia.fasta'
    multi_contigs = './flye_assembly/flye_multi_contig_loci.tsv'
    directory = './flye_assembly'
    
    ref_list = []
    
    for record in SeqIO.parse(reference, 'fasta'):
        ref_list.append([record.id, str(record.seq)])    

    flye_contig_list = pd.read_csv(multi_contigs, sep = '\t', index_col = 0)
    sample_list = list(flye_contig_list.columns.values)
    loci_list = list(flye_contig_list.index) 
    
    for sample in sample_list:
        
        print(sample)
        
        sample_contigs = pd.DataFrame(columns = ['LOCI', 'CONTIG', 'OUTCOME'])
        
        for loci in loci_list:
            if flye_contig_list[sample][loci]:
                    
                assemblyInfo = pd.read_csv(f'{directory}/{sample}/assembly/{loci}/flye/assembly_info.txt', sep = '\t', header = 0)
                
                contig_seqs = []
                
                for record in SeqIO.parse(f'{directory}/{sample}/assembly/{loci}/flye/assembly.fasta', 'fasta'):
                    contig_seqs.append([record.id, str(record.seq)])
                
                contig_seq = pd.DataFrame(contig_seqs, columns = ['CONTIG', 'SEQ'])
                
                temp_contigs = []

                for x in range(int(flye_contig_list[sample][loci])):
                    temp_contigs.append([loci, f'contig_{x + 1}', 'placeholder'])
                
                hit_table, qer_length = make_Alignment(directory, sample, loci, ref_list)
                
                try:
                
                    # hack to catch contigs that have multiple mappings
                    if hit_table['CONTIG'].nunique() == len(hit_table):
                    
                        for ctg in list(hit_table['CONTIG']):
                            index = int(ctg.split('_')[1]) - 1
                            temp_contigs[index][2] = 'mapped'
                        
                        contam_contigs = find_Containation(hit_table, list(assemblyInfo['#seq_name']), loci)
                        for ctg in contam_contigs:
                            index = int(ctg[0].split('_')[1]) - 1
                            temp_contigs[index][2] = ctg[1]
                        
                        partial_contigs =  find_Partial(hit_table, 0.75)
                        for ctg in partial_contigs:
                            index = int(ctg.split('_')[1]) - 1
                            temp_contigs[index][2] = 'partial'
                        
                        # not a fan of this solution but it will do...
                        remaining_contigs = [ctg for ctg in temp_contigs if ctg[2] == 'mapped']
                        remaining_contigs = [ctg[1] for ctg in remaining_contigs]
                        subset_hit_table = hit_table.loc[hit_table['CONTIG'].isin(remaining_contigs)]
                        
                        plot_Cigar(list(hit_table['CIGAR']), list(subset_hit_table['CONTIG']), qer_length, loci)
                        
                        if remaining_contigs:                                                                
                            
                            selected_hit = subset_hit_table.iloc[subset_hit_table['MAPQ'].idxmax()]
                            selected_idx = int(selected_hit['CONTIG'].split('_')[1]) - 1
        
                            temp_contigs[selected_idx][2] = 'selected'
                            
                            keep_seq = contig_seq.loc[contig_seq['CONTIG'] == selected_hit['CONTIG']]
                            
                            for contig in keep_seq.iterrows():
                                file = open(f'{directory}/{sample}/assembly/{loci}/flye/resolved.fasta', 'w')
                                file.write(f'>{loci}\n{contig[1]['SEQ']}\n')
                            
                        reject_seq = contig_seq.loc[contig_seq['CONTIG'] != selected_hit['CONTIG']]
                    
                        for contig in reject_seq.iterrows():
                            file = open(f'{directory}/{sample}/assembly/{loci}/flye/rejected.fasta', 'w')
                            file.write(f'>{contig[1]['CONTIG']}\n{contig[1]['SEQ']}\n')
                            
                except:
                    print(f'Skipped {loci}')
                
                temp_contigs = pd.DataFrame(temp_contigs, columns = ['LOCI', 'CONTIG', 'OUTCOME'])
                sample_contigs = pd.concat([sample_contigs, temp_contigs])
        
        sample_contigs.to_csv(f'{directory}/{sample}/{sample}_MAR.tsv', sep = '\t')
            
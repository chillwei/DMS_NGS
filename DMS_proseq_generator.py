#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 22:30:01 2024

@author: weiqiyao
"""
# this script contains the function to generate mock DMS library based on the wild type protein sequence including single amino acid replacement, insertion and deletion

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt
from Bio import pairwise2
import numpy as np
import os
from math import ceil
import itertools
from difflib import SequenceMatcher
import re
import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# function to generate single aa deletion DMS library
def single_aa_deletion(seq, start, end):
    """
    Perform a single amino acid deletion in the protein sequence.

    Parameters:
        seq (str): The wild type protein sequence.
        start position (int): The start position of the amino acid where we introduce deletion. (consider the first residue in the seq as 1)
        end position (int): The end position of the aa where we introduce deletion. 

    Returns:
        str: a list of mutated protein sequence.
    """
    mutseq_list = []
    for i in range(start , end+1,1):  
        mut_seq = seq[:(i-1)] + seq[i:]
        mutseq_list.append(mut_seq)
    #QC part to get rid of the duplicate mutants
    unique_mut  = list(set(mutseq_list))
    print('pre deduplicate library number is : ', len(mutseq_list))
    print('The diversity of single aa deletion library is: ', len(unique_mut))
    return unique_mut


# function to generate single aa inserrtion DMS library
def single_aa_insertion(seq, aa_list, start, end):
    """
    Perform a single amino acid insertion in the protein sequence.

    Parameters:
        aa_list: a list of amino acid you want to introduce as a insertion at each position for DMS
        seq (str): The wild type protein sequence.
        NOTE: insertion happens in between target aa position and target aa position -1
        start position (int): The start position of the amino acid where we introduce insertion. (consider the first residue in the protein sequence as 1)
        end position (int): The end position of the aa where we introduce insertion. 

    Returns:
        str: a list of mutated protein sequence.
    """
    mutseq_list = []
    for i in range(start , end+1,1):  
        for aa in aa_list : 
            mut_seq = seq[:(i-1)] + aa + seq[i-1:]
            mutseq_list.append(mut_seq)
    #QC part to get rid of the duplicate mutants
    unique_mut  = list(set(mutseq_list))
    print('pre deduplicate library number is : ', len(mutseq_list))
    print('The diversity of single aa insetrion library is: ', len(unique_mut))
    return unique_mut


# function to generate single aa replacement DMS library
def single_aa_replacement(seq, aa_list, start, end):
    """
    Perform a single amino acid replacement in each position of protein sequence except for the wild type residue. The region is defined by the input of start / end postion
    Parameters:
        aa_list: a list of amino acid you want to introduce as a replacement at each position for DMS
        seq (str): The wild type protein sequence.
        
        start position (int): The start position of the amino acid where we introduce replacement. (consider the first residue in the protein sequence as 1)
        end position (int): The end position of the aa where we introduce replacement. 

    Returns:
        str: a list of mutated protein sequence.
    """
    mutseq_list = []
    for i in range(start , end+1,1):  
        for aa in aa_list : 
            if aa != seq[i-1]: 
                mut_seq = seq[:(i-1)] + aa + seq[i:]
                mutseq_list.append(mut_seq)
    #QC part to get rid of the duplicate mutants
    unique_mut  = list(set(mutseq_list))
    print('pre deduplicate library number is : ', len(mutseq_list))
    print('The diversity of single aa replacement library is: ', len(unique_mut))
    return unique_mut
    

# profile the seq in NGS with the DMS designed seq to calculate the coverage of library diversity
# NOTE this function only finds whether the mutation exsit but not the exactly full protein sequence
def DMSlib_seq_profiling( NGS_proseq_list , DMS_proseq_list ):

    seq_not_capture = []
    #n = 0
    i = 0
    #Add 'HH' at the end of the seq
    #ins_smurfp_list_HHH = [string + 'HHH' for string in ins_smurfp_list]
    # find the extra seq in the dms oligos
    for proseqA in DMS_proseq_list:
        found = False
        # Iterate through listB
        for proseqB in NGS_proseq_list:
            # Check if stringA is a substring of stringB
            if proseqA in proseqB:
                found = True
                #n+=1
                #print(found,n)
                break
        # If stringA is not found in any stringB, append it to strings_only_in_listA
        if not found:
            i +=1
            seq_not_capture.append(proseqA)
            #print('not found', proseqA,i)
            
    coverage = (len(DMS_proseq_list)-i)/len(DMS_proseq_list) 
    
    # Print the strings that only exist in listA
    print("How many seqs that is not captured in DMS NGS lib :",i)
    print("coverage of the DMS is",coverage)

    return seq_not_capture


  # function to pull out the mutation information and convert it into a dataframe
def extract_mut_info(DMS_df, wtseq):
      df = DMS_df.copy()
      mut_loc_list = []
      mut_aa_list= []
      mut_type_list = []
      Proseq_list = []
      Proseqlen_list = []
      
      for idx, row in DMS_df.iterrows():
          mut_ID = row['ID']
          #mut_ID_list.append(mut_ID)
          DMS_seq = row['oligo_aa']
          mut_type = mut_ID.split('_')[1][:3]
          # annotate mutation type
          mut_type_list.append(mut_type)
           
          mut_loc = re.findall(r'\d+', mut_ID.split('_')[2])
          mut_loc_list.append(mut_loc[0])
          
          if mut_type == 'Del':
              del_aa = mut_ID.split('_')[2][:3]
              mut_aa_list.append(del_aa)
          
          else:
              mut_aa = mut_ID.split('_')[2][-3:]
              mut_aa_list.append(mut_aa)
          
          # recover full length of protein sequence for further DMS and NGS profiling
          mut_fragloc = re.findall(r'\d+', mut_ID.split('_')[1])
          mut_fragstart = int(mut_fragloc[0])
          mut_fragend = int(mut_fragloc[1])
          Proseq = wtseq[: (mut_fragstart-2)] + DMS_seq + wtseq[(mut_fragend +1) :]
          Proseq_list.append(Proseq)
          Proseqlen_list.append(len(Proseq))

    # Convert 3-letter abbreviations to single-letter abbreviations and join them
      #mut_aa_listv2 = [three_to_one[aa] for aa in mut_aa_list]
      df['Proseq'] = Proseq_list
      df['Proseq_len'] = Proseqlen_list
      df['mut_type'] = mut_type_list
      df['mut_aa'] = mut_aa_list
      df['mut_loc'] = mut_loc_list
      

      return df
         
         
         
# profile the seq in NGS with the DMS designed seq to calculate the coverage of library diversity
def NGS_DMS_capturedseq( NGS_proseq_list , DMS_proseq_list ):

    seq_capture = []
    #n = 0
    i = 0
    #Add 'HH' at the end of the seq

    # find the extra seq in the dms oligos
    for proseqA in DMS_proseq_list:
        #found = False
        # Iterate through listB
        for proseqB in NGS_proseq_list:
            # Check if stringA is a substring of stringB
            if proseqA in proseqB:
         #       found = True
                i +=1
                seq_capture.append(proseqB)
                #n+=1
                #print(found,n)
                #break
        # If stringA is not found in any stringB, append it to strings_only_in_listA
        
            #print('not found', proseqA,i)
            
    
    
    # Print the strings that only exist in listA
    print("How many seqs that is  captured in DMS NGS lib :",i)

    return seq_capture     


def _convert_DMSoligo2aa(df):
    Geneaa_list = []
    Geneaa_df = df.copy()
    for idx, row in df.iterrows():
        DNAseq = row['seq']
        DNAseqBsaI = DNAseq.replace('GGTCTC','@').replace('GAGACC','@')
        if DNAseqBsaI.count('@') != 2:
            print('Multiple BsaI found. Check the sequence manually ')
            print('ID:' + row['ID'])
            print('Sequence' + DNAseq)
        Genechunk = DNAseqBsaI.split('@')[1][2:-2]
        Gene_aa = Seq(Genechunk).translate()
        Geneaa_list.append(str(Gene_aa))
    

    Geneaa_df['oligo_aa'] = Geneaa_list
    return Geneaa_df
         


# dms_aa: Define all possible amino acids (adjust as needed)
# pass all possible amino acid designed in dms library and the dataframe generated above which represents counts of each unique variant seqinto the following function. Note: df has to be grouped by len for analysis.

def _missing_aa_heatmap(df, aa_list,mut_type, proseqlen):

    # Initialize a dictionary to hold the frequency data
    df2 = df[df['mut_type'] == mut_type]
    aa_matrix = {aa: [0]*proseqlen for aa in aa_list}
    # Calculate the frequency of each amino acid at each position
    for index, row in df2.iterrows():
        miss_aa = row['mut_aa']
        miss_loc = row['mut_loc']

        aa_matrix[miss_aa][int(miss_loc)] = 1


    # Convert the frequency matrix to a DataFrame for easier plotting
    aa_matrix_df = pd.DataFrame(aa_matrix)

    transposed_aa_matrix_df = aa_matrix_df.transpose()
    return transposed_aa_matrix_df


def _aa_distribution_heatmap(df, aa_list,mut_type, proseqlen):

    # Initialize a dictionary to hold the frequency data
    df2 = df[df['mut_type'] == mut_type]
    aa_matrix = {aa: [0]* (proseqlen + 1) for aa in aa_list}
    # Calculate the frequency of each amino acid at each position
    for index, row in df2.iterrows():
        mutaa = row['mut_aa']
        mutloc = row['mut_loc']
        mut_count = row['count']
        

        aa_matrix[mutaa][int(mutloc)] = mut_count


    # Convert the frequency matrix to a DataFrame for easier plotting
    aa_matrix_df = pd.DataFrame(aa_matrix).drop(index= 0)

    transposed_aa_matrix_df = aa_matrix_df.transpose()
    return transposed_aa_matrix_df



        
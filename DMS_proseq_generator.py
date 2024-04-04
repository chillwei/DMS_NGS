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
def extract_mut_info(DMS_df):
      df = DMS_df.copy()
      mut_loc_list = []
      mut_aa_list= []
      mut_type_list = []

      
      for idx, row in DMS_df.iterrows():
          mut_ID = row['ID']
          #mut_ID_list.append(mut_ID)
          DMS_seq = row['oligo_aa']
          mut_type = mut_ID.split('_')[1][:3]
         # annotate mutation type
          mut_type_list.append(mut_type)
          
          mut_loc = re.findall(r'\d+', mut_ID.split('_')[2])
          mut_aa = mut_ID.split('_')[2][-3:]
          mut_loc_list.append(mut_loc[0])
          mut_aa_list.append(mut_aa)

    # Convert 3-letter abbreviations to single-letter abbreviations and join them
      #mut_aa_listv2 = [three_to_one[aa] for aa in mut_aa_list]
      df['mut_type'] = mut_type_list
      df['mut_aa'] = mut_aa_list
      df['mut_loc'] = mut_loc_list

      return df
         
         
         
         
         
         
         











        
################
# Modules import
import pandas as pd
import numpy as np
import json
import os


################################################################################$$##########
# Reading information from json file. Used to extract the parameters from the `config.json`.
def read_json(path:str = "config.json") -> dict:
    """
    path : str -> path of the json file
    """

    with open('config.json') as config:
        config_f = json.load(config)

    return config_f


##########################################################################
# Creating folder according to the `congif.json` paths and program scheme.
def create_folders():
    config_dict = read_json()
    main_folders = [config_dict["input_folder"], config_dict["output_folder"]]
    subf_output = [config_dict["output_folder"] + i for i in ["//motif_figure", "//motif_data"]]

    for folder in main_folders + subf_output:
        if os.path.exists(folder) is False:
            os.mkdir(folder)
            print(f"> folder `{folder}` was created.")
        else:
            print(f"> folder `{folder}` exists, continuing.")


#########################################################
# Codons dictionary used for the nt sequence translation.
# Unified S" and S' codon -> for further processing 
codon_dic_updated = {
                'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                'TCT': "S", 'TCC': "S", 'TCA': "S", 'TCG': "S", # S' variant
                'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', # * for STOP
                'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',

                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', # S" variant

                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
                }


#########################################################
# Translating NT sequense and creating amino acid string.
def nt_transalte_104(dna_seq:str, 
                     aa_start:int = 1, 
                     aa_end:int = None) -> list:
    
    """
    cdr_seq: str -> nucleotide dna sequence format.
    aa_start: int -> first amino acid position.
    aa_end: int -> last amino acid position position.
    """
    
    # Validation of dna aligment and length (divided by 3).
    n_spaces = dna_seq.count("-")
    if  (n_spaces % 3 != 0) :
        raise Exception("Number of spacers dosent divide by 3, cheek sequence")
    
    # Getting dna seqs that divides by 3, removing right `edge` excess if needed.
    dna_seq = dna_seq[:len(dna_seq) - len(dna_seq) % 3]
    
    # Starting translation
    aa_translated = []
    aa_length = int((len(dna_seq)-len(dna_seq)%3)/3)

    # If defined specific aa positions
    if (aa_end != None):
        aa_length = aa_end
        pass

    # Translating the NT sequence
    for i in range(aa_start, aa_length):
        codon = dna_seq[i*3-3:i*3]
      
        if codon in list(codon_dic_updated.keys()):
            aa = codon_dic_updated[codon]  
        else:
            aa = "-"
         
        aa_translated.append(aa)
    
    # Returning string of amino acids from start to end
    return "".join(aa_translated)


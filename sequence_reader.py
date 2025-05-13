import pandas as pd
from time import time
import numpy as np

def read_fasta(file_path):
    df = pd.DataFrame()
    with open(file_path) as file:
        seq = ""
        for i, line in enumerate(file):                 
            if line.startswith(">"): # new entry starts
                if i > 0: # not the first entry in file
                    entry["SEQUENCE"] = seq # save last entries sequence
                    df = pd.concat([df, pd.DataFrame([entry])], ignore_index = True) # save last entry
                entry = {"SEQ_NAME": line[1:].rstrip()} # create new entry
                seq = ""
            else: seq += line.rstrip() # sequences are sometimes stretched over multiple lines

    entry["SEQUENCE"] = seq # save last entries sequence
    df = pd.concat([df, pd.DataFrame([entry])], ignore_index = True) # save last entry
    return df

def _ensembl_fasta_to_dict(fasta_entry):
    '''internal function to encode one line of fasta string into a dictionary'''
    _colons = np.unique([i for i, char in enumerate(fasta_entry) if char == ":"])
    stop = fasta_entry.find("description:")
    colons = _colons[_colons < stop]
    delim = np.unique([fasta_entry.rfind(" ",0, i) for i in colons])
    return {"SEQ_NAME":fasta_entry[0 :fasta_entry.find(" ")]} | {fasta_entry[delim[i]+1:colons[colons > delim[i]].min()] :  fasta_entry[colons[colons > delim[i]].min()+1:delim[i+1]] for  i in range(len(delim)-1)}


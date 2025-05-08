import pandas as pd
from time import time

def _ensembl_fasta_to_dict(fasta_entry, debug = False):
    '''internal function to encode one line of fasta string into a dictionary'''
    data = {}
    m1 = fasta_entry.find(" ")
    m2 = fasta_entry.find(" chromosome:",         m1)
    m3 = fasta_entry.find(" gene:",               m2)
    m4 = fasta_entry.find(" gene_biotype:",       m3)
    m5 = fasta_entry.find(" transcript_biotype:", m4)
    m6 = fasta_entry.find(" description:",        m5)
    m7 = fasta_entry.find(" [Source:",            m6)
    m8 = fasta_entry.find("]",                    m7)

    if debug: print(f"m1:{m1}, m2:{m2}, m3:{m3}, m4:{m4}, m5:{m5}, m6:{m6}, m7:{m7}, m8:{m8}", end = "\n\n")
    
    data["SEQ_NAME"]             = fasta_entry[   1 :m1]
    data["SEQ_TYPE"]             = fasta_entry[m1+1 :m2]
    data['CHROMOSOME']           = fasta_entry[m2+12:m3]
    data['GENE_ID']              = fasta_entry[m3+6 :m4]
    data['GENE_BIOTYPE']         = fasta_entry[m4+14:m5]
    data['TRANSCRIPT_BIOTYPE']   = fasta_entry[m5+20:m6]
    data['DESCRIPTION']          = fasta_entry[m6+13:m7]
    data['SEQUENCE']             = fasta_entry[m7+1:]
    
    return data


def read_ensembl_fasta(file_path, timeit=False):
    '''reads a fasta file from in the format of ensembl.org and returns as pandas.DataFrane'''
    if timeit: s = time()
    df = pd.DataFrame(columns = ['SEQ_NAME', 'SEQ_TYPE', 'CHROMOSOME', 'GENE_ID', 'GENE_BIOTYPE', 'TRANSCRIPT_BIOTYPE', 'DESCRIPTION', 'SEQUENCE'])
    buffer = ""
    with open(file_path) as file:
        buffer = ""
        for line in file:
            buffer += line.rstrip()
            z = buffer.rfind(">")
            if z>0:
                entry = buffer[:z].replace("\n", "")
                buffer = buffer[z:]
                df = pd.concat([df, pd.DataFrame([_ensembl_fasta_to_dict(entry)])], ignore_index=True)
                
    if timeit: print(f"{time()-s:.2g}s")
    return df


def read_microsynth_fasta(file_path):
    """reads the name and sequence from fasta file for a sample sent to microsynth seqlab gmbh for sequencing, returns as pandas.DataSeries"""
    with open(file_path) as file: content = [line for line in file]
    return pd.Series({"NAME": content[0][1:].rstrip(), "SEQUENCE": content[1]})



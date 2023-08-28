from collections import defaultdict
import os
import numpy as np


def get_mapped_res(uniprotid:str, PDB:str, chain:np.ndarray, PDB_to_UniProt:bool):

    """
    Creates a nested dictionary in the format of: 

    mapping[chain_id][UniProt_id] = [PDB_resid] (Parse UniProt resid into PDB resid)

    or 

    mapping[chain_id][PDB_resid] = [UniProt_id] (Parse PDB resid into UniProt resid)


    Format of annotation in mapping files:
    1. (Mapped) PDB_resid UniProt_resid
    2. (Insertion) PDB_resid Left_Empty
        a. Insertion is defined as gaps in subject seq but not in query seq
    3. (Missense) PDB_resid UniProt_resid
        a. Missense mutation is a type of Point mutation that changes 
        the identity of amino acid(s)
    4. (Deletion) UniProt_resid Left_Empty
    5. (Modified Residue) PDB_resid UniProt_resid
    denoted in the query sequence as 'X'

    """
    mapping_func = lambda: defaultdict(mapping_func)
    mapping = mapping_func()

    for chain_id in chain:

        if f"{uniprotid}-{PDB}_{chain_id}_resid_mapping.txt" in os.listdir("mapping/"):

            resid_mapping = np.loadtxt(f"mapping/{uniprotid}-{PDB}_{chain_id}_resid_mapping.txt", dtype=str)

            for i in range(resid_mapping.shape[0]):

                if not resid_mapping[i, 0].startswith("(Deletion)") and not resid_mapping[i, 0].startswith("(Insertion)"):

                    if PDB_to_UniProt: 

                        mapping[chain_id][int(resid_mapping[i, 1])] = [int(resid_mapping[i, 2])]

                    else: 
                        #UniProt_to_PDB
                        mapping[chain_id][int(resid_mapping[i, 2])] = [int(resid_mapping[i, 1])]
    return mapping

if __name__ == "__main__":

    # First unzip Ecoli_mapping.zip

    # Example
    print(get_mapped_res("P00448", "3OT7", ["B"], True)) # PDB to UniProt

    # print(get_mapped_res("P00448", "3OT7", ["B"], False)) # UniProt to PDB
   
   

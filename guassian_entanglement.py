import copy
import itertools
import os
import random
import typing
from collections import defaultdict
from operator import itemgetter
from warnings import filterwarnings
import numpy as np
from geom_median.numpy import compute_geometric_median  # used pip
from MDAnalysis import Universe
from numba import njit
from scipy.spatial.distance import cdist, pdist, squareform
from topoly import lasso_type  # used pip
import re

filterwarnings("ignore")

def pre_processing_pdb(pdb_file: str) -> None:

    """
    Pre-processing the PDB files (.pdb, .pdb1 (model 1) and alphafold pdbs) by removing everything after the last TER. 

    Function removes all model instances except the first one in the case of .pdb1 files. 
    
    Why remove the last TER?
        Residues after the TER are non-protein

    """

    pdb_data = np.loadtxt(f"PDBs/{pdb_file}", dtype=str, delimiter="\n")

    if pdb_file in os.listdir("Before_TER_PDBs/"):
        print(f"\033[4m{pdb_file}\033[0m is already processed. Please look in Before_TER_PDBs/")
        return

    if pdb_data[1].split()[1] == "ALPHAFOLD":
        print(f"ALPHAFOLD PDB \033[4m{pdb_file}\033[0m DOES NOT REQUIRE PROCESSING; COPYING PDB FILE TO Before_TER_PDBs/")
        os.system(f"cp PDBs/{pdb_file} Before_TER_PDBs/{pdb_file}")
        return
    
    else:
        print(f"PROCESSING \033[4m{pdb_file}\033[0m")

    last_TER_index = 0
    Model_indices = []

    for i in range(len(pdb_data)):

        if pdb_data[i].startswith("MODEL"):
            Model_indices.append(i)

    if len(Model_indices) > 1:

        for j in range(Model_indices[1]):

            if pdb_data[j].startswith("TER"):

                if j > last_TER_index:
                    last_TER_index = j

    elif len(Model_indices) == 1 or len(Model_indices) == 0:

        for j in range(len(pdb_data)):

            if pdb_data[j].startswith("TER"):

                if j > last_TER_index:
                    last_TER_index = j

    if last_TER_index == 0 and len(Model_indices) == 1:

        last_TER_index = len(pdb_data)
    
    elif last_TER_index == 0 and len(Model_indices) > 1:

        last_TER_index = Model_indices[1]

    with open(f"Before_TER_PDBs/{pdb_file}", "w") as f:

        for k in range(last_TER_index):
        
            f.write(f"{pdb_data[k]}\n")

@njit(fastmath=True)
def helper_dot(Runit: np.ndarray, dR_cross: np.ndarray) -> list:

    """
    Numba function to speed up dot product calculation. Ability 
    to use current GPU (if available) and CPU
    
    """

    return [np.dot(x,y) for x,y in zip(Runit,dR_cross)]

def point_60_rounding(num: float) -> float:

    """
    This function perform rounding to the nearest 0.6. 

    Ex: 0.61 rounds up to 1 but 0.59 rounds down to 0

    """
    if len(str(num).split("e")) == 2: 
        num = 0.0

    if num % 1 >= 0.60:
        rounded_num = round(num)
    else:
        rounded_num = int(str(num).split(".")[0])
    
    return rounded_num

def get_entanglements(coor: np.ndarray, l: int, termini_threshold: list, chain: str, resids: np.ndarray, resnames: np.ndarray) -> dict:

    """
    Find proteins containing non-covalent lasso entanglements.

    Entanglements are composed of loops (defined by native contacts) and crossing residue(s).

    """
    Nterm_thresh = termini_threshold[0]
    Cterm_thresh = termini_threshold[1]

    dist_matrix = squareform(pdist(coor))
    native_cmap = np.where(dist_matrix <= 8.0, 1, 0) 
    native_cmap = np.triu(native_cmap, k=4) 

    nc_indexs = np.stack(np.nonzero(native_cmap)).T
    
    range_l = np.arange(0, l-1)
    range_next_l = np.arange(1,l)

    coor = coor.astype(np.float32)
    R = 0.5*(coor[range_l] + coor[range_next_l])
    dR = coor[range_next_l] - coor[range_l]

    pair_array = np.asarray(list(itertools.product(dR,dR))) 

    x = pair_array[:,0,:]
    y = pair_array[:,1,:]

    dR_cross = np.cross(x, y)

    pair_array = np.asarray(list(itertools.product(R,R)))
    diff = pair_array[:,0,:] - pair_array[:,1,:]
    diff = diff.astype(np.float32)

    Runit = diff / np.linalg.norm(diff, axis=1)[:,None]**3 
    Runit = Runit.astype(np.float32)

    dot_matrix = helper_dot(Runit, dR_cross)
    dot_matrix = np.asarray(dot_matrix)
    dot_matrix = dot_matrix.reshape((l-1,l-1))

    nc_gdict = {} 

    for i,j in nc_indexs:

        loop_range = np.arange(i,j)
        nterm_range = np.arange(Nterm_thresh,i-5)
        cterm_range = np.arange(j+6,l-(Cterm_thresh + 1))

        gn_pairs_array = np.fromiter(itertools.chain(*itertools.product(nterm_range, loop_range)), int).reshape(-1, 2)
        gc_pairs_array = np.fromiter(itertools.chain(*itertools.product(loop_range, cterm_range)), int).reshape(-1, 2)

        if gn_pairs_array.size != 0:
            
            gn_vals = dot_matrix[gn_pairs_array[:,0],gn_pairs_array[:,1]]
            gn_vals = gn_vals[~np.isnan(gn_vals)] 
            gn_val = np.sum(gn_vals) / (4.0 * np.pi)
        
        else:
            gn_val = 0

        if gc_pairs_array.size != 0:
            
            gc_vals = dot_matrix[gc_pairs_array[:,0],gc_pairs_array[:,1]]
            gc_vals = gc_vals[~np.isnan(gc_vals)] 
            gc_val = np.sum(gc_vals) / (4.0 * np.pi)
        
        else:
            gc_val = 0
        
        rounded_gc_val = point_60_rounding(np.float64(abs(gc_val)))
        rounded_gn_val = point_60_rounding(np.float64(abs(gn_val)))

        if np.abs(rounded_gn_val) >= 1 or np.abs(rounded_gc_val) >= 1:
            
            nc_gdict[ (int(i), int(j)) ] = (gn_val, gc_val, rounded_gn_val, rounded_gc_val)

    missing_residues = find_missing_residues(resids)

    filtered_nc_gdict = loop_filter(nc_gdict, resids, missing_residues)

    entangled_res = find_crossing(coor.tolist(), filtered_nc_gdict, resids)

    filtered_entangled_res = crossing_filter(entangled_res, missing_residues)

    for ent in filtered_entangled_res:

        native_i = ent[0]
        native_j = ent[1]

        nc_i_resname = resnames[np.where(resids == native_i)][0]

        nc_j_resname = resnames[np.where(resids == native_j)][0]

        if nc_i_resname == "CYS" and nc_j_resname == "CYS":

            return None, []

    return filtered_entangled_res, missing_residues


def find_missing_residues(resids:np.ndarray) -> np.ndarray:

    """
    Find missing residues in pdb file

    """

    check_all_resids = np.arange(resids[0], resids[-1] + 1)

    missing_residues = np.setdiff1d(check_all_resids, resids)

    return missing_residues

def loop_filter(native_contacts: dict, resids: np.ndarray, missing_res: np.ndarray) -> dict:

    """
    Remove loops if there are three or more consecutive missing residues
    or the amount of any missing residues exceed 5% of the loop length 

    """

    for ij, values in native_contacts.items():

        native_i = resids[ij[0]]

        native_j = resids[ij[1]]

        rounded_gn = values[-2]

        rounded_gc = values[-1]

        check_loop = np.arange(native_i , native_j + 1) 

        loop_length = check_loop.size

        missing_res_loop = np.intersect1d(check_loop, missing_res)

        for index, diff_resid_index in itertools.groupby(enumerate(missing_res_loop), lambda ix : ix[0] - ix[1]):

            conseuctive_missing_residues = list(map(itemgetter(1), diff_resid_index))

            if len(conseuctive_missing_residues) >= 3 or len(missing_res_loop) > 0.05 * loop_length:

                native_contacts[ij] = None

    return native_contacts

def find_crossing(coor: np.ndarray, nc_data: dict, resids: np.ndarray) -> dict:

    """
    Use Topoly to find crossing(s) based on partial linking number

    """

    entangled_res = {}

    native_contacts = [[ij[0], ij[1]] for ij, values in nc_data.items() if values is not None]

    # reduction:
        # 1. each crossing must be 10 residues apart [default]
        # 2. first crossing should be at least 6 residues from the loop
        # 3. first crossing should be at least 5 residues from the closest termini

    data = lasso_type(coor, loop_indices=native_contacts, more_info=True, precision=0, density=0, min_dist=[10, 6, 5])
    # high precision, low denisty

    for native_contact in native_contacts:

        crossings = []

        native_contact = tuple(native_contact)

        if abs(nc_data[native_contact][-2]) >= 1: # if rounded_gn >= 1

            crossingN = [f"{cr[0]}{resids[int(cr[1:])]}" for cr in data[native_contact]["crossingsN"]]

            crossings += crossingN

        if abs(nc_data[native_contact][-1]) >= 1: # if rounded_gc >= 1
            
            crossingC = [f"{cr[0]}{resids[int(cr[1:])]}" for cr in data[native_contact]["crossingsC"]]
            
            crossings += crossingC

        gn = nc_data[native_contact][0]
        
        gc = nc_data[native_contact][1]

        ij_gN_gC = (resids[native_contact[0]], resids[native_contact[1]]) + (gn, gc) 

        entangled_res[ij_gN_gC] = np.unique(crossings)
        
    return entangled_res

def crossing_filter(entanglements: dict, missing_res: np.ndarray) -> dict:

    """
    Remove entanglements if there are any missing residues plus-and-minus 10 of the crossing(s)

    """
    
    for ij_gN_gC, crossings in entanglements.items():

        if crossings.size:

            check_crossing = []

            for crossing in crossings:

                reg_exp = re.split("\\+|-", crossing, maxsplit=1) # split the chirality

                check_crossing.append(np.arange(int(reg_exp[1]) - 10 , int(reg_exp[1]) + 11))

            check_crossing = np.concatenate(check_crossing)

            missing_res_cr = np.intersect1d(check_crossing, missing_res)

            if missing_res_cr.size:

                entanglements[ij_gN_gC] = None

    filtered_entanglements = {nc: re_cr for nc, re_cr in entanglements.items() if re_cr is not None and len(re_cr) > 0} 

    return filtered_entanglements

def calculate_native_entanglements(pdb_file: str) -> None:

    """
    Driver function that outputs native lasso-like self entanglements and missing residues for pdb and all of its chains if any

    """

    pdb = pdb_file.split(".")[0]

    if f"{pdb}_GE.txt" in os.listdir("unmapped_GE/"):
        print(f"\033[4m{pdb}\033[0m has non-covalent lasso entanglements. Please look in unmapped_GE/")
        return

    native_structure = Universe(f"Before_TER_PDBs/{pdb_file}", format="PDB")

    print(f"COMPUTING ENTANGLEMENTS FOR \033[4m{pdb}\033[0m")

    native_structure_dups = native_structure.select_atoms("name CA")

    termini_threshold = [5, 5]

    chains_to_analyze = set(native_structure_dups.segments.segids)

    for chain in chains_to_analyze:

        pdb_resids, i = np.unique(native_structure_dups.select_atoms(f"segid {chain}").resids, return_index=True)

        unique_native_structure = native_structure_dups.select_atoms(f"segid {chain}")[i]

        coor = unique_native_structure.positions

        chain_res = pdb_resids.size

        resnames = unique_native_structure.resnames

        if chain_res == 1:

            print(f"Skipping over chain {chain} for \033[4m{pdb}\033[0m since chain has only one alpha carbon")
            return

        if pdb_resids.size:

            ent_result, missing_residues = get_entanglements(coor, chain_res, termini_threshold, chain, pdb_resids, resnames)
            
            if ent_result: 

                for ij_gN_gC, crossings in ent_result.items():
                    if crossings.size:
                        with open(f"unmapped_GE/{pdb}_GE.txt", "a") as f:
                            f.write(f"Chain {chain} | ({ij_gN_gC[0]}, {ij_gN_gC[1]}, {crossings}) | {ij_gN_gC[2]} | {ij_gN_gC[3]}\n")

            if len(missing_residues):

                with open(f"unmapped_missing_residues/{pdb}_M.txt", "a") as f:
                    f.write(f"Chain {chain}: ")
                    for m_res in missing_residues:
                        f.write(f"{m_res} ")
                    f.write("\n")


if __name__ == "__main__":

    import multiprocessing as mp

    cores = len(os.sched_getaffinity(0))

    result_obj = set()

    directories = {"unmapped_GE", "unmapped_missing_residues", "Before_TER_PDBs", "clustered_unmapped_GE"}

    folder_exists = set(os.listdir("."))

    for folder in directories:
        if folder not in folder_exists:
            os.mkdir(f"{os.getcwd()}/{folder}") 

    for pdbs in os.listdir("PDBs/"):

        pre_processing_pdb(pdbs) 

    input_data = os.listdir("PDBs/")

    with mp.get_context("spawn").Pool(cores) as p:

        chunk = len(input_data) // (cores ** 2) + 1

        results = p.map_async(calculate_native_entanglements, iterable=input_data, chunksize=chunk)
        result_obj.add(results)
        
        for result in result_obj:
            result.get()

  

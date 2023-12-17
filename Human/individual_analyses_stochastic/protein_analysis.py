import os
import json
import numpy as np
import random
from collections import defaultdict
from MDAnalysis import Universe
from scipy.spatial.distance import cdist, pdist
import warnings
import multiprocessing as mp
from scipy.stats import rankdata
import ast
import sys
import pickle

warnings.filterwarnings("ignore")
                            
def null_dist_p_val_calc(uniprot: str):

    """
    Generate null distribution and p_value

    """
    pdb_resids = np.load("../DATA/pdb_resids_Human_V3.npz", allow_pickle = True)["arr_0"].tolist()
    pdb_coor = np.load("../DATA/pdb_coor_Human_V3.npz", allow_pickle = True)["arr_0"].tolist()
    functional_info = np.load("../DATA/human_V2_functional_protein.npz", allow_pickle = True)["arr_0"].tolist()
    entanglement_info = np.load("../DATA/human_non_covalent_lassos_4_5_no_knots.npz", allow_pickle = True)["arr_0"].tolist()

    print(uniprot)

    gene_pdb = list(pdb_resids[uniprot].keys())[0]

    rep_chain = list(functional_info[uniprot][gene_pdb].keys())[0]

    enriched_count = {}
    enriched_count[uniprot] = 0

    depleted_count = {}
    depleted_count[uniprot] = 0

    entangled_resids = []
    
    # 1. get the functional resids and entangled resids
    functional_resids = np.array(list(functional_info[uniprot][gene_pdb][rep_chain]), int)

    for rep_ent in entanglement_info[uniprot][gene_pdb][rep_chain]:

        entangled_resids.append(entanglement_info[uniprot][gene_pdb][rep_chain][rep_ent]["def_3"])
        # def_3 is the key for crossings

    entangled_resids = np.concatenate(entangled_resids)
    # remove np.unique() to keep observed state and random state the same

    neither_resids = set(pdb_resids[uniprot][gene_pdb][rep_chain]) - set(functional_resids) - set(entangled_resids)

    # 2. Grab the indices for those functional resids
    functional_indices = np.searchsorted(pdb_resids[uniprot][gene_pdb][rep_chain], functional_resids)
    entangled_indices = np.searchsorted(pdb_resids[uniprot][gene_pdb][rep_chain], entangled_resids)
    neither_indices = np.searchsorted(pdb_resids[uniprot][gene_pdb][rep_chain], np.array(list(neither_resids)))

    # 3. Perform the pairwise distances between every resid to every functional resdiues in protein 
    # Minimum is obtained by taking the local minima of each row in the below matrix
    # The matrix shape is (size of resids in protein, size of functional resids)
    # After taking the minimum the shape is (size of resids in protein, )

    minimum_dist = cdist(pdb_coor[uniprot][gene_pdb][rep_chain], pdb_coor[uniprot][gene_pdb][rep_chain][functional_indices]).min(axis=1) 
    
    # populate an array 
    data_labels = np.array([None] * len(pdb_resids[uniprot][gene_pdb][rep_chain]))

    # Use the populated array above and the indices to label accordingly
    data_labels[functional_indices] = "Functional" # functional residues
    data_labels[entangled_indices] = "Entangled" # entangled residues which can be subset of functional residues
    data_labels[neither_indices] = "Not_Entangled" # neither 
    
    combined_res_rank_labels = np.stack((pdb_resids[uniprot][gene_pdb][rep_chain], rankdata(minimum_dist), data_labels), axis=1)

    # Combine the resids, ranks and label like the following

    # [ [2, 1, "Entangled"], 
    #   [3, 2, "Entangled"], 
    #   [4, 3, "Entangled"], 
    #   [5, 4, "Entangled"], 
    #    .
    #    .
    #    . 
    #  [110, 109, 'Not_Entangled'] ]

    # get the corresponding inner arrays
    entangled_ranks = combined_res_rank_labels[entangled_indices]
    
    # sum up the ranks 
    observed_rank_sum = np.sum(entangled_ranks[:, 1])

    null_dist = []

    # population of all distance matrices
    nested_distances = np.load("../DATA/Array_of_distance_matrices.npz", allow_pickle = True)["arr_0"].tolist()

    gene_pdb_chain = f"{uniprot}_{gene_pdb}_{rep_chain}"

    with open("../DATA/human_sasa_dict_new.pkl", "rb") as reader:

        sasa_data_arrays = pickle.load(reader)[uniprot][gene_pdb][rep_chain]

    for iteration in range(50_000):

        ent_res = random_shuffling_w_spatial_correlation(pdb_resids, pdb_coor, entanglement_info, gene_pdb_chain, nested_distances, sasa_data_arrays)
        
        if ent_res is None:
            break
        
        rep_ents = list(ent_res[f"{uniprot}_{gene_pdb}_{rep_chain}"].keys())

        random_crossings = np.concatenate([ent_res[f"{uniprot}_{gene_pdb}_{rep_chain}"][rep_ent] for rep_ent in rep_ents])
        
        shuffled_combined_arr = combined_res_rank_labels[np.concatenate([np.where(r_cr == combined_res_rank_labels[:, 0])[0] for r_cr in random_crossings])]

        null_rank_sum = np.sum(shuffled_combined_arr[:, 1])

        if len(entangled_ranks) != len(shuffled_combined_arr) != len(random_crossings) != len(entangled_resids):

            print("something's wrong")
            sys.exit(0)

        null_dist.append(null_rank_sum)

        if null_rank_sum <= observed_rank_sum: 

            enriched_count[uniprot] += 1

        elif null_rank_sum >= observed_rank_sum:  

            depleted_count[uniprot] += 1  

    pval = {}
    pval[uniprot] = {}
    pval[uniprot]["Enrichment"] = 0
    pval[uniprot]["Depletion"] = 0
    pval[uniprot]["two-tailed"] = 0

    if ent_res is not None:

        null_dist.append(observed_rank_sum)

        pval[uniprot]["Enrichment"] = (enriched_count[uniprot] + 1) / 50_001
        pval[uniprot]["Depletion"] = (depleted_count[uniprot] + 1) / 50_001
        pval[uniprot]["two-tailed"] = min(pval[uniprot]["Enrichment"], pval[uniprot]["Depletion"]) * 2

        directory = "../DATA/individual_p_vals"

        func_site = "protein"
        func_site = func_site.replace(" ", "_")
        func_site = func_site.replace("-", "_")

        with open(f"{directory}/{func_site}/{uniprot}-{func_site}-Pval.json", "w") as js:
            js.write(json.dumps(pval))

    return pval

def check_random_position(difference_residues: set, gene_pdb_chain: str, sasa_data_arrays: dict):

    """
    Testing the random placement residue. 

    The random placement residue will give a probability of solvent exposure 

    if a random roll from [0, 1] including floats is > the probability of solvent exposre
    then pick a new random residue 

    otherwise select the random placement 

    """

    check_difference_residues = difference_residues.copy() # deep copy

    gene, rep_pdb, rep_chain = gene_pdb_chain.split("_")

    selected_placement = []

    while True:

        random_roll = np.random.rand()
    
        random_placement = np.random.choice(list(check_difference_residues)) 

        random_placement_idx = np.where(sasa_data_arrays["resids"] == random_placement)

        solvent_exposure = sasa_data_arrays["all_resids_avg_ratios"][random_placement_idx][0]

        random_placement_solv_exposure_probability = 0.0001 + 0.9358*np.exp(-256.3507*solvent_exposure) + 0.0519*np.exp(-22.7533*solvent_exposure) + 0.0123*np.exp(-2.6719*solvent_exposure)
        # bin = 0.01

        check_difference_residues.discard(random_placement)

        if random_roll <= random_placement_solv_exposure_probability:

            selected_placement.append(random_placement)
            break

        if len(check_difference_residues) == 0:

            # this means you went exhausted all the residues thus there no residues left to check/test
            break

    return selected_placement 

def random_shuffling_w_spatial_correlation(pdb_resids: dict, pdb_coor: dict, entanglement_info: dict, gene_pdb_chain: str, nested_dist: dict, sasa_data_arrays: dict):

    """
    Null Distribution Construction

    Criteria are:

        Maintains spatial orientation and distance

        Mimics the general location of the observed crossings 
        using a stochastic algorithm  
    
    """

    gene, gene_pdb, rep_chain = gene_pdb_chain.split("_")
    
    ent_res = {}

    ent_res[gene_pdb_chain] = {}

    for rep_ent in entanglement_info[gene][gene_pdb][rep_chain]:

        current_crossing_resids = [] 
        # resets for each entanglement. This means "with replacement" between entanglements

        ent_res[gene_pdb_chain][rep_ent] = []

        crossings = entanglement_info[gene][gene_pdb][rep_chain][rep_ent]["def_3"]

        num_of_cr = len(crossings)

        if num_of_cr > 1:

            random_dist_matrix = random.sample(nested_dist[num_of_cr], k = 1)[0]
            orginal_random_dist_matrix = random_dist_matrix.copy() # deep copy
            keep_track_of_dist_matrices = 1
        
        else:

            keep_track_of_dist_matrices = 0

        crossing_label = 0

        keep_track_first_res = set()

        keep_track_distances = set()

        while crossing_label < num_of_cr:

            if crossing_label == 0:

                if len(keep_track_first_res) == len(pdb_resids[gene][gene_pdb][rep_chain]) and len(keep_track_distances) == len(random_dist_matrix):

                    # pick a new distance matrix because exhuasted all first residue positons in protein
                    # and used all distances in the current matrix

                    random_dist_matrix = random.sample(nested_dist[num_of_cr], k = 1)[0]
                    orginal_random_dist_matrix = random_dist_matrix.copy() # deep copy

                    keep_track_of_dist_matrices += 1

                    keep_track_first_res.clear()
                    keep_track_distances.clear()

                difference = set(pdb_resids[gene][gene_pdb][rep_chain]) - keep_track_first_res # without replacement

                selected_placement = check_random_position(difference, gene_pdb_chain, sasa_data_arrays)

                if selected_placement:

                    keep_track_first_res.add(selected_placement[0]) 

                    ent_res[gene_pdb_chain][rep_ent].append(selected_placement[0]) 

                    current_crossing_resids.append(selected_placement[0])

                    crossing_label += 1

            else:

                distances = np.random.choice(random_dist_matrix, size = crossing_label, replace = False)

                keep_track_distances = keep_track_distances | set(distances)

                distances_i = np.concatenate([np.where(dist == random_dist_matrix)[0] for dist in distances])

                Thresholds = [[dist - 4, dist + 4] for dist in distances]

                potential_new_crossing = []

                # creating spheres from population of distance matrices
                for j, threshold in enumerate(Thresholds):

                    spherical_i = np.where(current_crossing_resids[j] == pdb_resids[gene][gene_pdb][rep_chain])[0]
                    
                    spherical_coor = pdb_coor[gene][gene_pdb][rep_chain][spherical_i][0]

                    spherical_dist = cdist([spherical_coor], pdb_coor[gene][gene_pdb][rep_chain])

                    layered_i = np.where((spherical_dist != 0) & (spherical_dist >= threshold[0]) & (spherical_dist <= threshold[1]))[1]

                    layered_res = pdb_resids[gene][gene_pdb][rep_chain][layered_i]

                    if layered_res.size:
                                
                        potential_new_crossing.append(layered_res)
                        # list of list

                    else: 
                        # if no such spheres that meet distance criteria;  
                        # pick a new first residue and restore the distance matrix

                        random_dist_matrix = orginal_random_dist_matrix

                        placed_cr_current = ent_res[gene_pdb_chain][rep_ent]

                        ent_res[gene_pdb_chain][rep_ent] = []

                        potential_new_crossing = []

                        crossing_label = 0

                        del current_crossing_resids[-len(placed_cr_current):]

                        break

                # if there are spheres that meet distance criteria ; start picking
                if potential_new_crossing:

                    inter = set(potential_new_crossing[0]).intersection(*potential_new_crossing)

                    difference = inter - set(current_crossing_resids) # without replacement

                    if inter and difference:

                        selected_placement = check_random_position(difference, gene_pdb_chain, sasa_data_arrays)

                        if selected_placement:

                            ent_res[gene_pdb_chain][rep_ent].append(selected_placement[0])

                            current_crossing_resids.append(selected_placement[0])

                            # delete the distance from distance matrix
                            random_dist_matrix = np.delete(random_dist_matrix, distances_i)

                            crossing_label += 1
                        
                    if not inter or not difference or not selected_placement:

                        # condition fails then
                        # pick a new first residue with the same distance matrix

                        random_dist_matrix = orginal_random_dist_matrix

                        placed_cr_current = ent_res[gene_pdb_chain][rep_ent]

                        ent_res[gene_pdb_chain][rep_ent] = []

                        crossing_label = 0

                        del current_crossing_resids[-len(placed_cr_current):]
                    
            if num_of_cr > 1 and keep_track_of_dist_matrices == len(nested_dist[num_of_cr]):
                # you exhausted all the distance matrices for that crossing

                print(f"exhausted all distance matrices for {rep_ent} and {gene_pdb_chain}")
                return

    return ent_res   
    
def main():
    
    cores = len(os.sched_getaffinity(0))

    result_obj = set()

    functional_sites = ["protein"]

    for func_site in functional_sites:

        print(func_site)

        entanglement_info = np.load("../DATA/human_non_covalent_lassos_4_5_no_knots.npz", allow_pickle = True)["arr_0"].tolist()
        
        functional_info = np.load("../DATA/human_V2_functional_protein.npz", allow_pickle = True)["arr_0"].tolist()

        if functional_info and entanglement_info:

            common_genes = set(functional_info.keys()).intersection(entanglement_info.keys())

            saved_func_site = func_site.replace(" ", "_")
            saved_func_site = saved_func_site.replace("-", "_")

            done_genes = set([pval_file.split("-")[0] for pval_file in os.listdir(f"../DATA/individual_p_vals/{saved_func_site}")])

            skipped = {"P16615", "Q8NB78", "O60551", "P30419", "Q9H0W9", "Q14562", "P11413", 
            "Q16790", "P14735", "O43570", "Q96HY7", "P04040", "Q9BZE1", "Q9NP92", "Q2NL82", 
            "Q9NR71", "Q8N1Q1", "P43166", "P27338", "Q9BU02", "P35219", 
            "Q96PN6", "P21399", "Q9UGM6", "Q06278", "P21397", "P16278", "Q92542"}

            common_genes = common_genes - done_genes - skipped

            print(len(common_genes))

            pdb_resids = np.load("../DATA/pdb_resids_Human_V3.npz", allow_pickle = True)["arr_0"].tolist()
            
            pdb_coor = np.load("../DATA/pdb_coor_Human_V3.npz", allow_pickle = True)["arr_0"].tolist()

            with mp.get_context("spawn").Pool(cores) as p:

                results = p.map_async(null_dist_p_val_calc, iterable = list(common_genes), chunksize = len(list(common_genes)) // (cores ** 2) + 1)
                
                pval  = results.get()
            

if __name__ == "__main__":

    # directories = {"../DATA/individual_p_vals", 
    #                 "../DATA/individual_p_vals/PDB_DNA_binding", 
    #                 "../DATA/individual_p_vals/PDB_free_ligand", 
    #                 "../DATA/individual_p_vals/PDB_metal_interface",
    #                 "../DATA/individual_p_vals/protein", 
    #                 "../DATA/individual_p_vals/PDB_RNA_binding", 
    #                 "../DATA/individual_p_vals/zinc_finger_region", 
    #                 "../DATA/individual_p_vals/active_site"}
    
    # folder_exists = set(os.listdir("."))
    # for folder in directories:
    #     os.makedirs(f"{os.getcwd()}/{folder}", exist_ok = True)

    main()


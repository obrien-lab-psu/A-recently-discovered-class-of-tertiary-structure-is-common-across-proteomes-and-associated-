import numpy as np
from collections import defaultdict
from scipy.spatial.distance import pdist
import ast
import sys

"""
Create distance matrices for the structural analysis

"""

all_entanglments = np.load("DATA/yeast_non_covalent_lassos_4_5_no_knots.npz", allow_pickle = True)["arr_0"].tolist()

pdb_resids = np.load("DATA/pdb_resids_Yeast_V3.npz", allow_pickle = True)["arr_0"].tolist()

pdb_coors = np.load("DATA/pdb_coor_Yeast_V3.npz", allow_pickle = True)["arr_0"].tolist()

QC = {}

single_cr = set()

mc_cr = set()

for each_gene in all_entanglments:

    num_cr_to_ents = defaultdict(list)

    rep_pdb = list(all_entanglments[each_gene].keys())[0]

    rep_chain = list(all_entanglments[each_gene][rep_pdb].keys())[0]

    gene_resid = pdb_resids[each_gene][rep_pdb][rep_chain]

    gene_coor = pdb_coors[each_gene][rep_pdb][rep_chain]

    count = 0

    for rep_ent in all_entanglments[each_gene][rep_pdb][rep_chain]:

        number_of_crossing = len(all_entanglments[each_gene][rep_pdb][rep_chain][rep_ent]["def_3"])

        # for a check
        count += number_of_crossing

        if number_of_crossing != 1:

            num_cr_to_ents[number_of_crossing].append(rep_ent)

    # check
    if count == len(all_entanglments[each_gene][rep_pdb][rep_chain]):
        single_cr.add((each_gene, rep_pdb, rep_chain))

    if num_cr_to_ents:

        for X_num_cr, ents in num_cr_to_ents.items():

            for ent in ents:

                crossings = ent[-X_num_cr:]

                crossings_coors = [gene_coor[np.where(cr == gene_resid)][0] for cr in crossings]

                # checked that len(crossings) == len(crossings_coors)

                if len(crossings) == len(set(crossings)):

                    pw_dist = pdist(crossings_coors)

                    if X_num_cr not in QC:
                        QC[X_num_cr] = {}

                    if (each_gene, rep_pdb, rep_chain) not in QC[X_num_cr]:

                        QC[X_num_cr][(each_gene, rep_pdb, rep_chain)] = defaultdict(list)
                    
                    QC[X_num_cr][(each_gene, rep_pdb, rep_chain)][ent].append(pw_dist)

X_num_cr_distances = defaultdict(list)

for num_of_cr, protein_w_cr_dist in QC.items():

    proteins = list(QC[num_of_cr].keys())

    for protein in proteins:

        # for a check
        mc_cr.add(protein)

        for ent, dist in protein_w_cr_dist[protein].items():

            X_num_cr_distances[num_of_cr].append(dist[0])

            # check
            n = num_of_cr

            expected_distances = (n**2 - n) / 2

            if expected_distances != len(dist[0]):

                print(protein, ent, dist, cr_check, len(dist[0]))

# np.savez("DATA/Array_of_distance_matrices.npz", X_num_cr_distances)







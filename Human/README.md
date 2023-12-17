
The purpose of this second README is to explain the code files and files in DATA/ for each organism folder

Run all Python files as python [name].py

structural_analysis.py --> answers the question are non-covalent entanglements (mainly crossing residues) enriched/depleted/neither near functional sites

- Requires human_non_covalent_lassos_4_5_no_knots.npz 
  - binary file containing a dictionary of all unique entanglements per representative gene 
- Requires human_V2_functional_all.npz
  - binary file containing a dictionary of all functional residues per representative gene for all functional sites
- Requires pdb_coor_Human_V3.npz
  - binary file containing a dictionary of all alpha-coordinates of representative structure for each gene 
- Requires pdb_resids_Human_V3.npz
  - binary file containing a dictionary of all pdb resids of representative structure for each gene
- Requires Array_of_distance_matrices.npz
  - binary file containing a dictionary of pairwise euclidean distances between crossing residues for all entanglements
- Requires human_sasa_dict_new.pkl
  - pickle file contains the probability distribution for the ratio of the surface accessible solvent area (SASA) of the crossings to the average SASA for the protein was generated with a bin size of 0.01 for each gene individually in the proteome
- **To open binary file: np.load("file.npz", allow_pickle=True)["arr_0"].tolist()**

individual_analysis.py --> answers the question are non-covalent entanglements enriched/depleted/neither near individual functional sites

- Requires human_non_covalent_lassos_4_5_no_knots.npz 
  - binary file containing a dictionary of all unique entanglements per representative gene 
- Requires individual functional sites binary files starting with "human_V2_functional_"
  - binary file containing a dictionary of all functional residues per representative gene for a specfic functional site
- Requires pdb_coor_Human_V3.npz
  - binary file containing a dictionary of all alpha-coordinates of representative structure for each gene 
- Requires pdb_resids_Human_V3.npz
  - binary file containing a dictionary of all pdb resids of representative structure for each gene
- Requires Array_of_distance_matrices.npz
  - binary file containing a dictionary of pairwise euclidean distances between crossing residues for all entanglements
- Requires human_sasa_dict_new.pkl
  - pickle file contains the probability distribution for the ratio of the surface accessible solvent area (SASA) of the crossings to the average SASA for the protein was generated with a bin size of 0.01 for each gene individually in the proteome

pval_correction_and_freq.py --> multiple hypothesis correction of p-values and display of structural anlaysis final result
  - Final result is final_human_corrected.png

create_pw_distance_matrices.py --> creates Array_of_distance_matrices.npz binary file

GO_enrichment_web.py --> Gene Ontology enrichment analysis and data collection code file

Human_UniProt_files_and_mapping/Human_Proteome_info.zip --> contains key information about each UniProt ID in Human Proteome such as PDB structures, functional residues from UniProt, function description, etc. 

Human_UniProt_files_and_mapping/Human_mapping.zip --> contains mapping files where the first column is annotation, second column is PDB numbering and last column is UniProt numbering
  - Please refer to parse_mapping.py for more description of the mapping files. This code parses these files and create a dictionary of either PDB to UniProt or UniProt to PDB 

DATA/Human_GO_stats_new.xlsx --> Output from GO_enrichment_web.py

DATA/rep_genes_with_entanglements.txt --> Representative genes with entanglements

DATA/percentages_of_knots_in_my_db.npz --> Proteins with Knots/slipknots. They are either knotted or not. 

DATA/representative_pdb_chain.txt --> Representative genes after filtering described in Methods

*Please note that all inputs thus far are in PDB numbering (entanglement and functional information)*

*Please note that all data obtained from web services such as UniProt and QuickGO (for Gene Ontology analysis) can change since they continuously update their databases.*

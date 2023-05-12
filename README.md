# Non-Covalent Lasso Entanglements

This repository contains code used in the following paper: doi. 

Specfically, it contains

    1. Code to generate non-covalent lasso entanglements given a PDB file (.pdb or .pdb1)
    2. Code and Data to run the Gene Ontology (GO) and Structural Enrichment statistical tests 
        a. They are located in the E.coli, Yeast and Humans folders. 
 
Please refer to the paper for the exact details of the Methods and Supplementary Information for the data. 

Note that the code was run and optimized on a linux machine.


## Run Locally

Clone the project

```bash
git clone https://github.com/VirajRana0/NCLE_GO_Project.git
```
Please make sure to install Python 3.8 or higher and the neccessary packages/modules using: 

```bash 
pip install -r requirements.txt
```

Optionally, you could either create a [Python virtual environment](https://docs.python.org/3/tutorial/venv.html) or a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands)

## Usage

To generate entanglements, please follow these steps: 

1. Put the protein files in the `PDBs` directory

2. Run the entanglement code

```python 
python guassian_entanglement.py
```
The output is located in `unmapped_GE` directory. 

Each file will have the same format where the columns represent:
- The chain
- Non-covalent entanglement in terms of (i, j, **r**) where **r** is a vector of crossing residues. Attached to each crossing residue is its chirality (either plus or negative).
- Linking number for the N terminus (GN)
- Linking number for the C terminus (GC)

`Example: Chain X | (i, j, [r₀, r₁, r₂, ..., rₙ]) | GN | GC`

*Please note that the computational time depends on the protein size and the number of CPU cores being used. The code automatically detects the numbers of CPU cores being used (only on linux).*

3. Cluster the entanglements and choose the correct cutoff for your species. Please see the code for the cutoffs. 

```python 
python clustering.py
```

The output is located in `clustered_unmapped_GE` and file format is similar as above except that vector **r** is part of the tuple (i, j, r₀, r₁, r₂, ..., rₙ)

*Please note that knots and slipknots are part of the code output. They should be filtered out as described in the paper.*

---

Code and Data for the enrichment tests are located in the Ecoli, Yeast and Human folders. Human individual_analysis is broken up into several scripts located in `Human/individual_analyses_stochastic` to reduce computational time.

```bash 
  GO_enrichment_web.py --> GO test
  structural_analysis.py --> Structural enrichment test
  individual_analysis.py ---> Structural enrichment test
```
## Entanglement Demo

An example has been calculated for the <u>E.coli 50S ribosomal subunit</u> (**PDB ID**: 6XZ7 chain S)

<p align="center">
  <img src="https://github.com/VirajRana0/NCLE_GO_Project/blob/main/img/6XZ7_S.gif" alt="My GIF">
</p>

Please see the VMD visualization state for 6XZ7 chain S under `img/` folder.

## Citation

If you found the code or data useful, please consider citing the paper: 

```bibtex 
@article {Non-covalent_lasso_entanglements,
  author       = {Viraj Rana, Ian Sitarik, Justin Petucci, Hyebin Song, and Edward P. O'Brien},
  journal      = {bioRxiv},
  title        = {A recently discovered class of tertiary structure is common across proteomes and associated with protein functions and biological processes},
  year         = {2023},
  doi          = {10.1101/2021.10.04.463034},
  URL          = {https://www.biorxiv.org/content/early/2021/10/04/2021.10.04.463034},
}
```

# Non-Covalent Lasso Entanglements

[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)

*Author: Viraj Rana*

Specfically, it contains

    1. Code to generate non-covalent lasso entanglements given a PDB file (.pdb or .pdb1)
    2. Code and Data to run the Gene Ontology (GO) and Structural Enrichment statistical tests 
        a. They are located in the E.coli, Yeast and Humans folders. 
 
Please refer to the paper for the exact details of the Methods and Supplementary Information for the data. 

Note that the code was run and optimized on a linux machine.


## Run Locally

Clone the project

```bash
git clone https://github.com/obrien-lab-psu/Non-covalent-Lasso-Entanglements-in-Folded-Proteins-Prevalence-Functional-Implications-and-Evolut.git
```
Please make sure to install Python version >= 3.9 and <=3.10. At the time of last update, the upper version limit was due to the Numba dependency. All the neccessary packages/modules can be installed using: 

```bash 
pip install -r requirements.txt
```

Optionally, you could either create a [Python virtual environment](https://docs.python.org/3/tutorial/venv.html) or a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands)

## Usage

To generate entanglements, please follow these steps under the `Entanglement_code` folder: 

1. Put the protein files in the `PDBs` directory

2. Run the entanglement code

```python 
python gaussian_entanglement_v4.5.py
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
  GO_enrichment_web.py --> Gene Ontology enrichment test
  structural_analysis.py --> Structural enrichment test
  individual_analysis.py ---> Structural enrichment test
  
  All the data required to run these tests are located under DATA/ folder for each organism. 
  They were obtained as described in the Methods.
```
### Note for the Structural enrichement tests. 
All the p-values have been calculated already, so running the code will not do anything.
Why? The code check the `DATA` folder if a p-value has been calculated. If so, the code will not run. 
Please provide new input if you plan to use the code. 

## Entanglement Demo

An example has been calculated for the <u>E.coli 50S ribosomal subunit</u> (**PDB ID**: 6XZ7 chain S)

<p align="center">
  <img src="https://github.com/obrien-lab-psu/Non-covalent-Lasso-Entanglements-in-Folded-Proteins-Prevalence-Functional-Implications-and-Evolut/blob/main/img/6XZ7_S.gif" alt="My GIF">
</p>

Please see the VMD visualization state for 6XZ7 chain S under `img/` folder. 

## Citation

If you found the code or data useful, please consider citing the paper: 

```bibtex 
@article {Non-covalent_lasso_entanglements,
  author       = {Viraj Rana, Ian Sitarik, Justin Petucci, Yang Jiang, Hyebin Song, and Edward P. O'Brien},
  journal      = {Journal of Molecular Biology},
  title        = {Non-covalent Lasso Entanglements in Folded Proteins: Prevalence, Functional Implications, and Evolutionary Significance},
  year         = {2023},
  doi          = {https://doi.org/10.1016/j.jmb.2024.168459},
  URL          = {https://www.sciencedirect.com/science/article/pii/S0022283624000251},
}
```

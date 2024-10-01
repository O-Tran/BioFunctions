# Bioinformatic Functions
### Overview
BioFunction is a Python repository that contains scripts for analyzing biomolecular structures. This repository includes two primary scripts:

1. ```sortInteraction.py```: Downloads PDB files, analyzes ligand-protein interactions and compiles the results into a CSV file.
2. ```averageBfactor.py```: Calculates the average B-factor (or PLDDT score) of a protein structure from a PDB file.

### Features
- Download PDB files directly from the RCSB database.
- Analyze various interactions, including hydrogen bonds, hydrophobic interactions, and pi-stacking.
- Calculate average B-factors or PLDDT scores from protein structures.
- Output results in a structured CSV format.

### Requirements
Ensure you have the following Python packages installed:

     pip install pandas biopython plip

## Scripts
#### 1. ```sortInteraction.py```- Analyze interactions between ligands and proteins in PDB structures.
#### Usage:
    python3 sortInteraction.py -l [list_of_PDBs] -o [output.csv] 
#### Example:
    python3 sortInteraction.py -l 4WNM 4OGS 6NBL -o output.csv 

#### 2. ```averageBfactor.py```- Calculate the average B-factor (or PLDDT score) of a given PDB structure.

#### Usage:
    python3 averageBfactor.py [input.pdb]

#### Example:
    python3 averageBfactor.py 1BNA.pdb

## License
This repository is open-source. Feel free to use and modify the scripts for your projects. Contributions are welcome!

## Getting Started
### Clone the repository:


    git clone https://github.com/O-Tran/BioFunction.git
    cd BioFunction

### Install the required packages:

    pip install pandas biopython plip

Run the desired script according to the provided usage instructions.

## Contribution
If you'd like to contribute to this repository, please create a pull request with your proposed changes, or open an issue for discussion.

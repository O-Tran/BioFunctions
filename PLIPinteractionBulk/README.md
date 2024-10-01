# Python script to "bulk-analyze" the interaction between ligands and proteins

## Overview
The sortInteraction.py script downloads PDB files, analyzes the interactions between ligands and proteins, and compiles the results into a CSV file. It extracts interaction details such as hydrogen bonds, hydrophobic interactions, and pi-stacking, saving these interactions in a structured format for further analysis.

### Version
- Version: 1.0.0
- Date: April 26, 2023
- Author: Oanh TN Tran
### Requirements
Python Packages
    - pandas: For handling data frames and CSV output.
    - urllib: To download PDB files and the SMILES database.
    - plip: Used to analyze protein-ligand interactions.


You can install the required packages using:

``` pip install pandas plip biopython ```

### Usage
---
To run the script, you need to provide a list of PDB codes and an output file name for the CSV results.

### Command-Line Usage
``` python3 sortInteraction.py -l [list_of_PDBs] -o [output.csv]  ```

### Examples
    
    # Example 1: Using multiple PDB codes

    python3 sortInteraction.py -l 4WNM 4OGS 6NBL -o output.csv 

    # Example 2: Providing PDB codes as strings

    python3 sortInteraction.py -l '4WNM' '4OGS' '6NBL' -o output.csv 

### Input:
- List of PDB IDs: A space-separated list of PDB codes that you want to analyze.
- Output file: The name of the output CSV file that will store the interaction results.
### Output:
- A CSV file containing ligand-protein interaction details, including hydrogen bonds, hydrophobic interactions, and pi-stacking interactions.
-----

### Function Descriptions
``` WritePDB(pdb_id) ```

Downloads the PDB file using its PDB ID and stores it in the local pdbfiles directory.

```ParsePDB(inputpdb)```

Parses the PDB file using the plip library and extracts the interaction information for each ligand. Returns a dictionary of interactions by site.

```downloadsmilesdata()```

Downloads the latest SMILES database from the RCSB website for ligand-to-SMILES conversion.

```ligPDB2smiles(ligand_PDB, inputdir)```

Converts a PDB ligand ID (3-letter code) to its corresponding SMILES string using the downloaded SMILES database.

```isligand(lig, inputdir)```

Determines if the given PDB ligand is a valid ligand based on its SMILES string.

```getresidue(item)```

Generates a string that identifies a specific residue based on atom and chain data.

```pdb2dictligand(pdbfile)```

Creates a dictionary mapping ligand atoms to their respective PDB file atom IDs.

```sortInteraction(pdb)```

Main function that combines all the interaction data and outputs it as a pandas DataFrame.

## Notes
- The script automatically downloads an updated SMILES database for ligand-to-SMILES conversion. This file can be deleted after use if not needed.

- If a PDB file contains no ligands, the script will notify the user and skip that file.
---
### License
Feel free to use and modify this script for your own projects.

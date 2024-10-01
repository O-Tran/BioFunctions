# PLDDT/B-Factor Score Calculator
## Overview
This Python script calculates the average PLDDT or B-factor score of a given protein structure from a PDB file. The PLDDT score is typically used for AlphaFold models, while the B-factor score is used to assess the atomic displacement in X-ray crystallography structures.

### Prerequisites
The script requires Python and the Biopython package to run.

### Python Packages:
Biopython: Used for parsing and processing PDB files.
Install Biopython with the following command:

``` pip install biopython ```


### How to Use the Script

You can run the script by passing a PDB file as an argument. The script will parse the structure, calculate the average B-factor (or PLDDT) for each residue, and output the overall average score for the entire structure.

### Usage:
```python script_name.py <input_pdb_file> ```

### Example:
``` python script_name.py example.pdb```


This will output the average B-factor or PLDDT score for the provided PDB structure.

### Input:
- A valid PDB file (e.g., a .pdb file that contains protein structure data).
### Output:
- The overall average B-factor or PLDDT score of the structure.

If the structure does not contain B-factor or PLDDT scores, the script will print a message notifying you that the PDB file lacks this information.

### Command-line Arguments:
- pdb: The path to the PDB file to analyze.

### Code Breakdown
1. findBfactor(structurefile):

    - Parses the input PDB file and extracts the B-factor (or PLDDT) scores from the atoms of each residue.
    - Averages the scores for each residue and calculates the overall average for the entire structure.
    - Prints and returns the average score.
2. Command-line Interface:

    - The script uses argparse to handle command-line input. 
    - You provide the path to the PDB file, and the script processes it.


## Output Example:
``` PLDDT / Bfactor of this structure is  -- 45.87 ```

### Error Handling:
If the PDB file lacks any B-factor or PLDDT scores, the script will notify you with the following message:

``` <filename> doesn't have Bfactor or PLDDT score ```


### License
Feel free to modify and use this script for your own projects.

from Bio.PDB.PDBParser import PDBParser
import argparse
import statistics

""" 
Run this script to gather PLDDT / bfactor scores total average of the structure 
- input: Structure PDB 
- output: Bfactor 

Python need: 
Biopython   


"""

def findBfactor(structurefile):  
    ## Parse pdb into biopython 
    structure = PDBParser(PERMISSIVE = False , QUIET=False).get_structure('X',structurefile)
    residue_list = structure.get_residues() ## its actually a generator , that get emptied when cycled
    
    
    ## go to individaul residues to get the individual atoms 
    totalBfactor = [] # Clear 
    for res in residue_list:
        ResAtom  = res.get_atoms()
        
        ## collect average Bfactor/PlDDT from these atoms 
        bfactorlist = []  # clear 
        for atom in ResAtom: 
            bfactorlist.append(atom.get_bfactor()) 
        
        ## Average out and add to the total residues number Bfactor/PlDDT 
        AveragePerRes= statistics.mean(bfactorlist) 
        bfactorlist = []  # Clear 
        totalBfactor.append(AveragePerRes) 
        
        
    # average all the residues's Bfactor/PlDDT 
    AverageProtBfactor = statistics.mean(totalBfactor) 
    totalBfactor = []  # clear 
    print ("PLDDT / Bfactor of this structure is  -- {}".format(AverageProtBfactor)) 
    return AverageProtBfactor



    
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Get PLDDT Scores from PDBS and/or Bfactor scores when you input the PDB")

    parser.add_argument('pdb', help='Input pdb file here')  ### or points towards directory file 
    args = parser.parse_args()
    avgbfactor = findBfactor(args.pdb)
    if avgbfactor == 0: 
        print("{args.pdb} doesn't have Bfactor or PLDDT score")


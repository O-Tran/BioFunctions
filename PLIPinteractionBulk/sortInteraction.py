import sys 
import os 
import pandas as pd
import urllib.request
from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport

"""
sortInteraction.py

Objective: take PDB code, download PDB files, look into the different interactions and compile them into pandas/CSV output 

Version: 1.0.0   
Date: April 26 2023 
Author: Oanh TN Tran 

use with harmonic environment - 
Usage: python3 sortInteraction.py -l [listpdb] -o [output.csv]
Example: python3 sortInteraction.py -l 4WNM 4OGS 6NBL -o output.csv
Example: python3 sortInteraction.py -l '4WNM' '4OGS' '6NBL' -o output.csv 

Note: This script also download an updated Smiles Database that can be convert to smiles (this can be deleted later) 
"""

def WritePDB (pdb_id):  #convert PDB id to files
    """
    Converts pdb number to local files 
    Input: PDB id only in str ('') 
    Return: name of the file in str format  
    """
    downloadurl='http://files.rcsb.org/download/'
    
    datadir = os.getcwd()+"/pdbfiles"
    
    if os.path.exists(datadir) == False:
        os.mkdir("./pdbfiles")
    
    pdbfn = pdb_id + ".pdb"
    print(pdbfn)
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)    
    
    try: 
        path=urllib.request.urlretrieve(url, outfnm) 
    #fetch then get it to that directory
        return outfnm
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None
    
def ParsePDB (inputpdb): #parsePDB complex into PLIP and get dictionary of interaction 
    """
    Add in PDB (directoryfile) and check for biding site 
    inputfile: str for location of pdbfile
    Return : Dict of dictionary of the sites and the individual interaction 
    """
    #parsing in the PDB complex
    protComplex = PDBComplex()
    protComplex.load_pdb(inputpdb) #works! 

    #determeine if there's ligand, if not, return None 
    if len(protComplex.ligands) == 0: 
        print (inputpdb + " - has no ligands")
        return None 
    
    #if have ligand, Characterize the different binding properties of ligand 
    #this includes, but not limited to: Hbond, hydrophobic interaction, and pistacking 
    else:
        for lig in protComplex.ligands: 
            protComplex.characterize_complex(lig)
        
        #make dictonary for the sites 
        sites = {} 
        for key, site in sorted(protComplex.interaction_sets.items()):
            binding_site = BindingSiteReport(site) 
            #print(binding_site.hbond_features)
            keys = (
                "hbond",
                "hydrophobic",
                "pistacking"
                )
            interactions = {k: getattr(binding_site, k + "_info") for k in keys} 
            sites [key] = interactions
        return sites #checked 

def downloadsmilesdata(): #download Smiles txt data 
    """ download database from RSBC online """
    smiles_url = 'http://ligand-expo.rcsb.org/dictionaries/Components-smiles-stereo-oe.smi'
    outputdir = "SmileDictionary.txt" 
    
    urllib.request.urlretrieve(smiles_url, outputdir)
    return outputdir 

def ligPDB2smiles(ligand_PDB, inputdir): #convert Ligand ID to smiles 
    """
    Converts pdb ligand 3 letter code to smiles via RDkit  
    Input: PDB id only in str ('') 
    Return: smiles in str format  
    """
    data= pd.read_csv(inputdir, sep="\t", header=None)
    data.columns = ['smiles', 'PDBcode', 'name']

    #data.to_dict('index') #works! yay 
    smilesdict = data.set_index('PDBcode').to_dict() #works now 
    smile_out = smilesdict['smiles'].get(ligand_PDB)

    return smile_out 

def isligand(lig, inputdir): #inputdir from smiles.txt downloaded in updated* downloadsmilesdata (input str, output str) 
    smile = ligPDB2smiles(lig, inputdir) 
    if smile.find("C") + smile.find("c") == 0: 
        print(lig + " is not a ligand" )
        return None
    else:
        return lig, smile

def getresidue (item): #get the amino acids in the interaction (input dict, output str)
    Residue = str(item[1])+str(item[0])
    return Residue

def pdb2dictligand(pdbfile): #create dictionary of atoms in ligand and what that atomtype corresponds to (input file str, output dict)
    ligatomdict = {} 
    with open (pdbfile, 'r') as pdbtxt: 
        read = pdbtxt.readlines()
        for line in read:
            splitline = line.split("\n'")
            doublesplit = splitline[0].split(' ')
            if doublesplit[0] == "HETATM":
                ligatomdict[int(doublesplit[1])] = doublesplit[3]
    return ligatomdict

def sortInteraction(pdb): #Combinging together! output forms a dataframe 
    smilesdata = downloadsmilesdata()
    maindf = pd.DataFrame()    
    pdbdir = WritePDB(pdb)
    if pdbdir == None:
        print('{} does not exist'.format(pdb))
        return None 
    else:
        sitedf = pd.DataFrame()
        sites = ParsePDB(pdbdir)
        if sites ==None:
            #print ('{} does not have ligand'.format(pdb))
            return None 
        for site in sites:
            
            acceptor = [] 
            donor = [] 
            hydrophobic = [] 
            pistacking = [] 

            lig, chain, resID = site.split(':')
            smiles = ligPDB2smiles(lig, smilesdata)
            ligatomdict = pdb2dictligand(pdbdir)

            if len(smiles) == 0:
                pass
            else: 
                #hydrogenbond interactions 
                hbonds = sites[site]['hbond']
                for i in range(len(hbonds)): 
                    residue = getresidue(hbonds[i])

                    #identify if Ligand is donor or acceptor 
                    if hbonds[i][10] == True: #protein is donor, ie ligand is acceptor 
                        Ligatom = ligatomdict[hbonds[i][13]] #14 
                        acceptor.append((residue, Ligatom) )
                    elif hbonds[i][10] == False: #protein is NOT donor, ie ligand is donor  
                        Ligatom = ligatomdict[hbonds[i][11]] #12 
                        donor.append((residue, Ligatom))

                #hydrophobic interactions 
                hydrophobics = sites[site]['hydrophobic']
                for i in range(len(hydrophobics)):
                    residue = getresidue(hydrophobics[i]) 
                    Ligatom = ligatomdict[hydrophobics[i][7]]
                    hydrophobic.append((residue, Ligatom))


                #pistacking interactions
                pistackings = sites[site]['pistacking']
                for i in range(len(pistackings)):
                    residue = getresidue(pistackings[i])
                    Ligatom = ligatomdict[pistackings[i][11]]
                    pistacking.append((residue, Ligatom))

                data = {'PDB': pdb , 'Smiles': smiles,
                       'Hbond Acceptor':[acceptor],
                       'Hbond Donor':[donor], 
                       'Hydrophobic':[hydrophobic],
                       'PiStacking':[pistacking]}

                df1 = pd.DataFrame(data, columns=data)
            maindf = pd.concat([maindf,df1])
            #maindf.to_csv(outputCSV) 
    return maindf

#'''
if __name__ == "__main__": #parsing argument into python code (Adjust to list later )
    import argparse
    parser = argparse.ArgumentParser() 
    parser.add_argument('-l','--list', nargs="+", help='<Input list of PDB> Set flag', required=True)
    parser.add_argument('-o','--output', help="Output name/directory")

    args = parser.parse_args()
    
    #sortInteraction(args.list, args.o)
    
    #run function with parsed argument 
    df =pd.DataFrame() 
    for pdb in args.list:
        outputdf = sortInteraction(pdb) 
        df=pd.concat([df, outputdf])
    df.to_csv(args.output)
    

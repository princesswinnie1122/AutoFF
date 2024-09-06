import re
import json
from rdkit import Chem
from rdkit.Chem import BRICS

def fragment_molecule(smiles):
    """Fragment a molecule using BRICS."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, ["Invalid SMILES string."]
    
    frags = BRICS.BreakBRICSBonds(mol)
    frag_smiles = [Chem.MolToSmiles(frag, True) for frag in Chem.GetMolFrags(frags, asMols=True)]
    return frag_smiles

def remove_num(smiles):
    """
    Replace BRICS attachment points like [16*], [14*], etc.
    with [*], removing the numbers.
    """
    # Replace [number*] patterns with [*]
    return re.sub(r'\[\d+\*\]', '[*]', smiles)

def search_database(fragments, database_file="BB_dic.json"):
    """Replace linkers, normalize fragments, and check if they exist in the database."""
    mod_frag = [remove_num(frag) for frag in fragments]
    frag_status = []

    with open(database_file, 'r') as db_file:
        database = json.load(db_file)
    
    for frag in mod_frag:
        found = False
        for bb_name, bb_smiles_list in database.items():
            if frag in bb_smiles_list:
                found = True
                frag_status.append((bb_name, frag, True))
                break
        if not found:
            frag_status.append((None, frag, False))
    
    return frag_status
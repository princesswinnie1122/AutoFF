import re
import json
from rdkit import Chem
from brics import BRICS
from clean import clean_num, clean_all

def fragment_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, ["Invalid SMILES string."]
    
    frags = BRICS.BreakBRICSBonds(mol)
    frag_smiles = [Chem.MolToSmiles(frag, True) for frag in Chem.GetMolFrags(frags, asMols=True)]
    return frag_smiles

def clean_num(smiles):
    return re.sub(r'\[\d+\*\]', '[*]', smiles)


def search_database(fragments, bb_dic_file="BB_dic.json", bb_info_file="BB_info.json"):
    frag_status = []

    with open(bb_dic_file, 'r') as db_file:
        bb_dic = json.load(db_file)
    
    with open(bb_info_file, 'r') as info_file:
        bb_info = json.load(info_file)

    for frag in fragments:
        frag = clean_num(frag) # [*]
        #print(frag)

        found_in_dic = False
        found_in_info = False
        bb_name = None

        mol_star = Chem.MolFromSmiles(frag)
        cans_frag_star = Chem.MolToSmiles(mol_star, isomericSmiles=True, canonical=True)
        cans_frag_star = re.sub(r'\*', '[*]', cans_frag_star)

        for name, bb_smiles_list in bb_dic.items():
            if cans_frag_star in bb_smiles_list:
                found_in_dic = True
                bb_name = name
                break

        if found_in_dic:
            frag_status.append((bb_name, cans_frag_star, True))
        else:
            #cleaned = clean_all(frag)
            cleaned = re.sub(r'\[\*\]', '[H]', frag)

            #print(cleaned)
            mol_h = Chem.MolFromSmiles(cleaned)
            cans_frag_h = Chem.MolToSmiles(mol_h, isomericSmiles=True, canonical=True)

            for name, info in bb_info.items():
                if clean_all(info['smiles']) == cans_frag_h:
                    found_in_info = True
                    bb_name = name
                    break

            if found_in_info:
                frag_status.append((bb_name, cans_frag_star, "Missing Connectors"))
            else:
                frag_status.append((None, cans_frag_star, "Missing SMILES"))

    return frag_status
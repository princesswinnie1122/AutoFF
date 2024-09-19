import re
import json
from rdkit import Chem
from rdkit.Chem import BRICS
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

'''
def search_database(fragments, bb_dic_file="BB_dic.json", bb_info_file="BB_info.json"):
    frag_status = []

    with open(bb_dic_file, 'r') as db_file:
        bb_dic = json.load(db_file)
    
    with open(bb_info_file, 'r') as info_file:
        bb_info = json.load(info_file)

    for frag in fragments:
        frag = clean_num(frag)
        print(frag)

        found_in_dic = False
        found_in_info = False
        bb_name = None

        for name, bb_smiles_list in bb_dic.items():
            if frag in bb_smiles_list:
                found_in_dic = True
                bb_name = name
                break

        if found_in_dic:
            frag_status.append((bb_name, frag, True))
        else:
            cleaned = clean_all(frag)
            print(cleaned)
            mol = Chem.MolFromSmiles(cleaned)
            Chem.Kekulize(mol)
            cans_frag = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True, canonical=True)
            #cans_frag = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            for name, info in bb_info.items():
                if clean_all(info['smiles']) == cans_frag:
                    found_in_info = True
                    bb_name = name
                    break

            if found_in_info:
                frag_status.append((bb_name, frag, "Missing Connectors"))
            else:
                frag_status.append((None, frag, "Missing SMILES"))

    return frag_status
'''

def search_database(fragments, bb_dic_file="BB_dic.json", bb_info_file="BB_info.json"):
    frag_status = []

    with open(bb_dic_file, 'r') as db_file:
        bb_dic = json.load(db_file)
    
    with open(bb_info_file, 'r') as info_file:
        bb_info = json.load(info_file)

    for frag in fragments:
        frag = clean_num(frag)
        print(f"Cleaned fragment: {frag}")

        found_in_dic = False
        found_in_info = False
        bb_name = None

        for name, bb_smiles_list in bb_dic.items():
            if frag in bb_smiles_list:
                found_in_dic = True
                bb_name = name
                break

        if found_in_dic:
            frag_status.append((bb_name, frag, True))
        else:
            cleaned = clean_all(frag)
            print(f"Fully cleaned fragment: {cleaned}")
            
            try:
                mol = Chem.MolFromSmiles(cleaned)
                if mol is None:
                    raise ValueError("Invalid SMILES")
                
                Chem.Kekulize(mol)
                cans_frag = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True, canonical=True)
            except Exception as e:
                print(f"Error processing SMILES: {e}")
                frag_status.append((None, frag, f"Invalid SMILES: {str(e)}"))
                continue

            for name, info in bb_info.items():
                if clean_all(info['smiles']) == cans_frag:
                    found_in_info = True
                    bb_name = name
                    break

            if found_in_info:
                frag_status.append((bb_name, frag, "Missing Connectors"))
            else:
                frag_status.append((None, frag, "Missing SMILES"))

    return frag_status
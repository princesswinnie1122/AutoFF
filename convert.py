import json
import re
from rdkit import Chem 
from rdkit.Chem import AllChem as rdkit
import stk
from clean import clean_para


def replace(smiles):
    return re.sub(r'Br|I|F', '[*]', smiles)

# To make sure that the user is entering the SMILES they want
def info_to_dic(input_file='BB_info.json', output_file='BB_dic.json'):
    with open(input_file, 'r') as f:
        data = json.load(f)

    BB_dic = {}
    problematic_bbs = []

    for molecule_name, molecule_info in data.items():
        for BB in data[molecule_name]["BB"]:
            try:
                molecule = molecule_name
                smiles = clean_para(data[molecule]["smiles"])
                BB_name = BB

                # Validate SMILES
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError(f"Invalid SMILES: {smiles}")

                Br_num = data[molecule]['BB'][BB_name].get('Br')
                I_num = data[molecule]['BB'][BB_name].get('I')
                F_num = data[molecule]['BB'][BB_name].get('F')

                building_block = stk.BuildingBlock(smiles)
                mol_rdkit =  building_block.to_rdkit_mol()
                atom = list(mol_rdkit.GetAtoms())

                for n, a in enumerate(atom):
                    if Br_num is not None and n+1 in Br_num:
                        a.SetAtomicNum(35)  # Br
                    if I_num is not None and n+1 in I_num:
                        a.SetAtomicNum(53)  # I 
                    if F_num is not None and n+1 in F_num:
                        a.SetAtomicNum(9)   # F

                smiles = Chem.MolToSmiles(mol_rdkit, True)
                canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles, True))
                modified_smiles = replace(canonical_smiles)
                BB_dic[BB_name] = [modified_smiles]  # Keep it as a list for consistency

            except Exception as e:
                print(f"Error processing {BB_name} ({molecule_name}): {str(e)}")
                problematic_bbs.append((BB_name, molecule_name, str(e)))
                continue

    with open(output_file, 'w') as json_file:
        json.dump(BB_dic, json_file, indent=4)

    if problematic_bbs:
        print("\nThe following building blocks could not be processed:")
        for bb, molecule, error in problematic_bbs:
            print(f"- {bb} ({molecule}): {error}")

    return BB_dic

if __name__ == "__main__":
    info_to_dic()
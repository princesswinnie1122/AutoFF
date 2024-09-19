import re

def clean_para(smiles):
    cleaned = re.sub(r'\(\)', '', smiles)
    cleaned = re.sub(r'=\)', ')', cleaned)
    cleaned = re.sub(r'\(=', '(', cleaned)
    return cleaned

def clean_all(smiles):
    cleaned = re.sub(r'\[\d*\*\]', '', smiles)
    cleaned = clean_para(cleaned)
    return cleaned

def clean_num(smiles):
    return re.sub(r'\[\d+\*\]', '[*]', smiles)
import re

def normalize_fragment(smiles):
    """Remove BRICS attachment points from a SMILES string."""
    return re.sub(r'\[\d+\*\]', '', smiles)

def check_fragments_in_database(fragments, database):
    """Check if normalized fragments exist in the database."""
    normalized_fragments = [normalize_fragment(frag) for frag in fragments]
    fragment_status = {}
    
    for frag in normalized_fragments:
        found = False
        for entry in database.values():
            if frag == entry["smiles"]:
                found = True
                break
        fragment_status[frag] = found
    
    return fragment_status

# Example usage
fragments = [
    "[14*]c1cc([14*])nc(C)n1",
    "[16*]c1cc([16*])cc([16*])c1",
    "[16*]c1cccnc1",
    "[16*]c1ccncc1"
]

database = {
    "Benzene": {
        "smiles": "c1ccccc1",
        "BB": {
            "1BA": {"Br": [7]},
            # More building blocks
        }
    },
    "Pyridine": {
        "smiles": "c1ccncc1",
        "BB": {
            "2BA": {"Br": [7]},
            # More building blocks
        }
    }
}

status = check_fragments_in_database(fragments, database)
print(status)

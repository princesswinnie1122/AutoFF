from openbabel import pybel

def get_bond_connectivity(smiles):
    """Generate molecular structure and print bond connectivity."""
    mol = pybel.readstring("smi", smiles)
    mol.addh()  # Add hydrogens
    mol.make3D()  # Generate 3D structure

    print(f"Atom connectivity for: {smiles}")
    for bond in mol.OBMol.Bonds():
        # Get atom indices (OpenBabel uses 1-based indexing)
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()
        bond_order = bond.GetBondOrder()

        atom1_symbol = mol.OBMol.GetAtom(atom1).GetType()
        atom2_symbol = mol.OBMol.GetAtom(atom2).GetType()

        print(f"{atom1_symbol} (atom {atom1}) - {atom2_symbol} (atom {atom2}), Bond order: {bond_order}")

# Example usage:
smiles = "Brc1ccccc1"  # Example: Bromobenzene
get_bond_connectivity(smiles)

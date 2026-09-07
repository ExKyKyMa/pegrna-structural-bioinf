from rdkit import Chem


smiles = "C(CC(=O)O)CN=[N+]=[N-]"
molecule = Chem.MolFromSmiles(smiles)

if molecule is None:
    print("Invalid SMILES")
else:
    print("Valid SMILES")
    print(Chem.MolToSmiles(molecule))
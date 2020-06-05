"""
Creates a png of a list of SMILES
Sorts the structures by taking first structure and then going from most similar to least similar
"""
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors

data = open('data/benzodiazepine_activator/total_smiles.txt', 'r')

molecules = []

for line in data.readlines():
    compound = Chem.MolFromSmiles(line[:-1])
    molecules.append((rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound), compound))

sorted_mols = sorted(molecules[1:], key=lambda x: DataStructs.DiceSimilarity(molecules[0][0], x[0]))

img = Draw.MolsToGridImage([x[1] for x in sorted_mols], molsPerRow=10, subImgSize=(200, 200))
img.save('molecules.png')
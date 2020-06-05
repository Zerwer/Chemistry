"""
Creates a png of a list of SMILES
Iterates over compounds, finds most similar then finds most similar to that compound until all compounds are arranged
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


def sort_similarity(mols, sort):
    largest = [0, None]
    for mol in mols:
        similarity = DataStructs.DiceSimilarity(sort[len(sort)-1][0], mol[0])
        if similarity > largest[0]:
            largest = [similarity, mol]

    mols.remove(largest[1])

    sort.append(largest[1])
    if len(mols) > 0:
        sort_similarity(mols, sort)
    return sort


sorted_mols = sort_similarity(molecules[1:], [molecules[0]])
img = Draw.MolsToGridImage([x[1] for x in sorted_mols], molsPerRow=10, subImgSize=(200, 200))
img.save('molecules.png')

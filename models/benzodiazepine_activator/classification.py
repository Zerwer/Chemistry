# Simple benzodiazepine classfication model based on similarity
#   to other benzodiazepines
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
from sklearn.model_selection import train_test_split

data = open('data/benzodiazepine_activator/total_smiles.txt', 'r')

molecules = []

for line in data.readlines():
    compound = Chem.MolFromSmiles(line[:-1])
    molecules.append((rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound), compound))


def model(mol, active_mols):
    similarities = []
    for active_mol in active_mols:
        similarities.append(DataStructs.DiceSimilarity(mol[0], active_mol[0]))

    return max(similarities)


train_mols, test_mols = train_test_split(molecules, test_size=0.1, random_state=1)

for test_mol in test_mols:
    prediction = model(test_mol, train_mols)
    print(prediction)

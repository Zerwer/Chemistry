"""
All general functions go here instead of chemical models
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs


# Function that returns n amount of similar molecules and their pKa
def pka_similarities(mol, mol_set, n):
    mol_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(mol))
    similarity = []
    for molecule in mol_set:
        sim = DataStructs.DiceSimilarity(mol_fp, molecule[2])
        similarity.append([sim, molecule[1]])

    return np.asarray(sorted(similarity)[:n]).flatten()

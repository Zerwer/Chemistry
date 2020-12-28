# Common reusable functions
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs, Descriptors


# Consumes a molecule and a list of molecules and produces a list of n most
#   similar molecules from mol_set
# Length of mol_set <= n
def pka_similarities(smile, mol_set, n):
    mol = Chem.MolFromSmiles(smile)
    mol_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
    similarity = []
    for molecule in mol_set:
        sim = DataStructs.DiceSimilarity(mol_fp, molecule[2])
        similarity.append([sim, molecule[1]])

    return np.asarray(sorted(similarity)[:n]).flatten()


# Converts the logS value of a molecule to mg/mL
def logs_to_mg_ml(logs, mol):
    mg_l = (10.0 ** float(logs)) * float(Descriptors.ExactMolWt(mol))
    return mg_l

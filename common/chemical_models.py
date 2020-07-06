"""
Removes redundancy from models that use other models as an input
"""
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors, DataStructs
import numpy as np


# All models inherit this as all the models are pickled
class ChemicalModel:
    def __init__(self, name):
        self.model = pickle.load(open('run_models/' + name + '_model.pkl', 'rb'))
        self.scaler = pickle.load(open('run_models/' + name + '_scaler.pkl', 'rb'))


class LogP(ChemicalModel):

    def run(self, fingerprint):
        return self.model.predict(self.scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]


class LogPSolubility(ChemicalModel):

    def run(self, logP):
        return self.model.predict(self.scaler.transform(np.asarray(logP).reshape(1, -1)))[0]


class AtomPairSolubility(ChemicalModel):

    def run(self, fingerprint):
        return self.model.predict(self.scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]


class MeltingPoint(ChemicalModel):

    def run(self, sol, logP):
        return self.model.predict(self.scaler.transform(np.asarray([sol, logP]).reshape(1, -1)))[0]


class AcidpKa(ChemicalModel):

    def run(self, fingerprint):
        return self.model.predict(self.scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]


# Function that returns n amount of similar molecules and their pKa
def pka_similarities(mol, mol_set, n):
    mol_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(mol))
    similarity = []
    for molecule in mol_set:
        sim = DataStructs.DiceSimilarity(mol_fp, molecule[2])
        similarity.append([sim, molecule[1]])

    return np.asarray(sorted(similarity)[:n]).flatten()


class AcidSimilarity(ChemicalModel):

    def run(self, mol, acid_set):
        sim = pka_similarities(mol, acid_set, 512)
        return self.model.predict(self.scaler.transform(np.asarray(sim).reshape(1, -1)))[0]


class CombinedSolubility(ChemicalModel):

    def run(self, mol, logP, logP_sol, atom_pair_sol):
        mw = Descriptors.ExactMolWt(mol)
        rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        ap = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[a]'))) / mol.GetNumHeavyAtoms()
        esol = 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * rb - 0.74 * ap

        return self.model.predict(self.scaler.transform(np.asarray([logP_sol, atom_pair_sol, esol]).reshape(1, -1)))[0]

# Removes redundancy from models that use other models as an input
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from functions import pka_similarities


# All models inherit this as all the models are pickled
class ChemicalModel:
    def __init__(self, name):
        self.model = pickle.load(open('run_models/' + name +
                                      '_model.pkl', 'rb'))
        self.scaler = pickle.load(open('run_models/' + name +
                                       '_scaler.pkl', 'rb'))


class LogP(ChemicalModel):

    def run(self, fingerprint):
        scaled = self.scaler.transform(np.asarray(fingerprint).reshape(1, -1))
        return self.model.predict(scaled)[0]


class LogPSolubility(ChemicalModel):

    def run(self, logP):
        scaled = self.scaler.transform(np.asarray(logP).reshape(1, -1))
        return self.model.predict(scaled)[0]


class AtomPairSolubility(ChemicalModel):

    def run(self, fingerprint):
        scaled = self.scaler.transform(np.asarray(fingerprint).reshape(1, -1))
        return self.model.predict(scaled)[0]


class MeltingPoint(ChemicalModel):

    def run(self, sol, logP):
        scaled = self.scaler.transform(np.asarray([sol, logP]).reshape(1, -1))
        return self.model.predict(scaled)[0]


class AcidpKa(ChemicalModel):

    def run(self, fingerprint):
        scaled = self.scaler.transform(np.asarray(fingerprint).reshape(1, -1))
        return self.model.predict(scaled)[0]


class AcidSimilarity(ChemicalModel):

    def run(self, mol, acid_set):
        similarity = pka_similarities(mol, acid_set, 512)
        scaled = self.scaler.transform(np.asarray(similarity).reshape(1, -1))
        return self.model.predict(scaled)[0]


class CombinedSolubility(ChemicalModel):

    def run(self, mol, logP, logP_sol, atom_pair_sol):
        mw = Descriptors.ExactMolWt(mol)
        rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        ap = (len(mol.GetSubstructMatches(Chem.MolFromSmarts('[a]')))
              / mol.GetNumHeavyAtoms())
        # Formula for estimating solubility (ESOL)
        esol = 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * rb - 0.74 * ap

        combined = np.asarray([logP_sol, atom_pair_sol, esol]).reshape(1, -1)
        scaled = self.scaler.transform(combined)
        return self.model.predict(scaled)[0]

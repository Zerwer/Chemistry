# Removes redundancy from models that use other models as an input
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors, MACCSkeys
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
from functions import pka_similarities, logs_to_mg_ml


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


class GeneralPKa(ChemicalModel):

    def run(self, fingerprint):
        scaled = self.scaler.transform(np.asarray(fingerprint).reshape(1, -1))
        return self.model.predict(scaled)[0]


class AcidpKa(GeneralPKa):
    pass


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


class AllProperties:
    def __init__(self, logP, sol, mp, pka):
        self.logP = logP
        self.sol = sol
        self.mp = mp
        self.pka = pka


class AllModels:
    descriptors = {'logP': 'logP',
                   'solubility': 'sol',
                   'melting': 'mp',
                   'pKa': 'pka'}

    def __init__(self, logP, logP_sol, atom_pair_sol,
                 combined_sol, melting_point, pKa):
        self.logP_model = LogP(logP)
        self.logP_solubility_model = LogPSolubility(logP_sol)
        self.atom_pair_sol_model = AtomPairSolubility(atom_pair_sol)
        self.combined_model = CombinedSolubility(combined_sol)
        self.melting_point_model = MeltingPoint(melting_point)
        self.pKa_model = GeneralPKa(pKa)

    def predict(self, mol):
        fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
        avalon = GetAvalonFP(mol)
        maacs = MACCSkeys.GenMACCSKeys(mol)

        logP = self.logP_model.run(fp)
        logP_sol = self.logP_solubility_model.run(logP)
        atom_pair_sol = self.atom_pair_sol_model.run(fp)
        combined_sol = self.combined_model.run(mol, logP,
                                                      logP_sol, atom_pair_sol)
        mg_ml_sol = logs_to_mg_ml(combined_sol, mol)
        mp = self.melting_point_model.run(combined_sol, logP)
        pka = self.pKa_model.run(avalon + maacs + fp)

        return AllProperties(logP, mg_ml_sol, mp, pka)

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
                   'pKa': 'pka',
                   'molecular weight': 'mol_wt'}

    def __init__(self, logP, logP_sol, atom_pair_sol,
                 combined_sol, melting_point, pKa):
        self.logP_model = LogP(logP)
        self.logP_solubility_model = LogPSolubility(logP_sol)
        self.atom_pair_sol_model = AtomPairSolubility(atom_pair_sol)
        self.combined_model = CombinedSolubility(combined_sol)
        self.melting_point_model = MeltingPoint(melting_point)
        self.pKa_model = GeneralPKa(pKa)

    # Function predicts molecules descriptors from selected_descriptors
    #   Since some descriptors have redundant calculations the list options
    #   is used to avoid redundancy
    def predict(self, mol, selected_descriptors):
        options = [0, 0, 0, 0, 0]
        return_properties = {}

        for option in selected_descriptors:
            if option == 'logP':
                options[0] = 1
            elif option == 'sol':
                options[0] = 1
                options[1] = 1
            elif option == 'mp':
                options[0] = 1
                options[1] = 1
                options[2] = 1
            elif option == 'pka':
                options[3] = 1
            elif option == 'mol_wt':
                options[4] = 1

        fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)

        if options[0]:
            logP = self.logP_model.run(fp)
            return_properties['logP'] = logP

        if options[1]:
            logP_sol = self.logP_solubility_model.run(logP)
            atom_pair_sol = self.atom_pair_sol_model.run(fp)
            combined_sol = self.combined_model.run(mol, logP,
                                                   logP_sol, atom_pair_sol)
            mg_ml_sol = logs_to_mg_ml(combined_sol, mol)
            return_properties['sol'] = mg_ml_sol

        if options[2]:
            mp = self.melting_point_model.run(combined_sol, logP)
            return_properties['mp'] = mp

        if options[3]:
            avalon = GetAvalonFP(mol)
            maacs = MACCSkeys.GenMACCSKeys(mol)
            pka = self.pKa_model.run(avalon + maacs + fp)
            return_properties['pka'] = pka

        if options[4]:
            wt = rdMolDescriptors.CalcExactMolWt(mol)
            return_properties['mol_wt'] = wt

        return return_properties

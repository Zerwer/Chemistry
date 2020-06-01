"""
Removes redundancy from models that use other models as an input
"""
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
import numpy as np
import os

print(os.getcwd())



def init(name):
    model = pickle.load(open('run_models/' + name + '_model.pkl', 'rb'))
    scaler = pickle.load(open('run_models/' + name + '_scaler.pkl', 'rb'))

    return [model, scaler]


def run_logP(model, fingerprint):
    return model[0].predict(model[1].transform(np.asarray(fingerprint).reshape(1, -1)))[0]


def run_logP_sol(model, logP):
    return model[0].predict(model[1].transform(np.asarray(logP).reshape(1, -1)))[0]


def run_atom_pair_solubility(model, fingerprint):
    return model[0].predict(model[1].transform(np.asarray(fingerprint).reshape(1, -1)))[0]


def run_combined_solubility(model, mol, logP, logP_sol, atom_pair_sol):
    mw = Descriptors.ExactMolWt(mol)
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    ap = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[a]'))) / mol.GetNumHeavyAtoms()
    esol = 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * rb - 0.74 * ap

    return model[0].predict(model[1].transform(np.asarray([logP_sol, atom_pair_sol, esol]).reshape(1, -1)))[0]
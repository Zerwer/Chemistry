# Predicts solubility using the combined solubility model
# Usage: python3 solubility.py SMILES
import pickle, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors

logP_model = pickle.load(open('run_models/logP_model.pkl', 'rb'))
logP_scaler = pickle.load(open('run_models/logP_scaler.pkl', 'rb'))

logP_solubility_model = pickle.load(open('run_models/logS_logP_model.pkl', 'rb'))
logP_solubility_scaler = pickle.load(open('run_models/logS_logP_scaler.pkl', 'rb'))

solubility_model = pickle.load(open('run_models/water_solubility_model.pkl', 'rb'))
solubility_scaler = pickle.load(open('run_models/water_solubility_scaler.pkl', 'rb'))

combined_model = pickle.load(open('run_models/combined_solubility_model.pkl', 'rb'))
combined_scaler = pickle.load(open('run_models/combined_solubility_scaler.pkl', 'rb'))

compound = Chem.MolFromSmiles(sys.argv[1])

fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound)

logP = logP_model.predict(logP_scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]

logP_sol = logP_solubility_model.predict(logP_solubility_scaler.transform(np.asarray(logP).reshape(1, -1)))[0]

sol = solubility_model.predict(solubility_scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]

mw = Descriptors.ExactMolWt(compound)
rb = rdMolDescriptors.CalcNumRotatableBonds(compound)
ap = len(compound.GetSubstructMatches(Chem.MolFromSmarts('[a]'))) / compound.GetNumHeavyAtoms()
esol = 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * rb - 0.74 * ap

combined = combined_model.predict(combined_scaler.transform(np.asarray([logP_sol, sol, esol]).reshape(1, -1)))[0]


def logs_conversion(logS):
    return (10.0**float(logS)) * float(Descriptors.ExactMolWt(Chem.MolFromSmiles(sys.argv[1])))


print('LogP solubility prediction: ' + str(logs_conversion(logP_sol)))
print('Atomic Pair solubility prediction: ' + str(logs_conversion(sol)))
print('ESOL solubility prediction: ' + str(logs_conversion(esol)))
print('Combined solubility prediction: ' + str(logs_conversion(combined)))
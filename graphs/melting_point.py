import pickle
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors

data = open('data/melting_point/mp.txt', 'r')

logP_model = pickle.load(open('run_models/logP_model.pkl', 'rb'))
logP_scaler = pickle.load(open('run_models/logP_scaler.pkl', 'rb'))

logP_solubility_model = pickle.load(open('run_models/logS_logP_model.pkl', 'rb'))
logP_solubility_scaler = pickle.load(open('run_models/logS_logP_scaler.pkl', 'rb'))

solubility_model = pickle.load(open('run_models/water_solubility_model.pkl', 'rb'))
solubility_scaler = pickle.load(open('run_models/water_solubility_scaler.pkl', 'rb'))

combined_model = pickle.load(open('run_models/combined_solubility_model.pkl', 'rb'))
combined_scaler = pickle.load(open('run_models/combined_solubility_scaler.pkl', 'rb'))

melting_point_model = pickle.load(open('run_models/melting_gse_model.pkl', 'rb'))
melting_point_scaler = pickle.load(open('run_models/melting_gse_scaler.pkl', 'rb'))

x = []
y = []

for line in data.readlines():
    split = line.split(' ')
    compound = Chem.MolFromSmiles(split[0])

    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound)

    logP = logP_model.predict(logP_scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]

    logP_sol = logP_solubility_model.predict(logP_solubility_scaler.transform(np.asarray(logP).reshape(1, -1)))[0]

    sol = solubility_model.predict(solubility_scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]

    mw = Descriptors.ExactMolWt(compound)
    rb = rdMolDescriptors.CalcNumRotatableBonds(compound)
    ap = len(compound.GetSubstructMatches(Chem.MolFromSmarts('[a]'))) / compound.GetNumHeavyAtoms()
    esol = 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * rb - 0.74 * ap

    combined = combined_model.predict(combined_scaler.transform(np.asarray([logP_sol, sol, esol]).reshape(1, -1)))[0]

    mp = melting_point_model.predict(melting_point_scaler.transform(np.asarray([combined, logP]).reshape(1, -1)))[0]

    x.append(float(mp))
    y.append(float(split[1][:-1]))

plt.scatter(x, y, s=1)
plt.xlabel('Predicted')
plt.ylabel('Experimental')
plt.title('Melting Point Prediction')
plt.show()

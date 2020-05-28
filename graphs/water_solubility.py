# Plot different neural networks to predict solubility
import pickle
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors

data = open('data/water_solubility/aqsol.txt', 'r')

logP_model = pickle.load(open('run_models/logP_model.pkl', 'rb'))
logP_scaler = pickle.load(open('run_models/logP_scaler.pkl', 'rb'))

logP_solubility_model = pickle.load(open('run_models/logS_logP_model.pkl', 'rb'))
logP_solubility_scaler = pickle.load(open('run_models/logS_logP_scaler.pkl', 'rb'))

solubility_model = pickle.load(open('run_models/water_solubility_model.pkl', 'rb'))
solubility_scaler = pickle.load(open('run_models/water_solubility_scaler.pkl', 'rb'))

combined_model = pickle.load(open('run_models/combined_solubility_model.pkl', 'rb'))
combined_scaler = pickle.load(open('run_models/combined_solubility_scaler.pkl', 'rb'))

x1, x2, x3, x4 = [], [], [], []
y1, y2, y3, y4 = [], [], [], []

for line in data.readlines():
    split = line.split(' ')

    compound = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound)
    logP = logP_model.predict(logP_scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]
    logP_sol = logP_solubility_model.predict(logP_solubility_scaler.transform(np.asarray(logP).reshape(1, -1)))[0]
    sol = solubility_model.predict(solubility_scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]

    # Additional ESOL empirical model to increase accuracy
    mw = Descriptors.ExactMolWt(compound)
    rb = rdMolDescriptors.CalcNumRotatableBonds(compound)
    ap = len(compound.GetSubstructMatches(Chem.MolFromSmarts('[a]')))/compound.GetNumHeavyAtoms()
    esol = 0.16 -  0.63*logP - 0.0062*mw + 0.066*rb - 0.74*ap

    combined = combined_model.predict(combined_scaler.transform(np.asarray([logP_sol, sol, esol]).reshape(1, -1)))[0]

    x1.append(float(logP_sol))
    x2.append(float(sol))
    x3.append(float(esol))
    x4.append(float(combined))

    y1.append(float(split[1][:-1]))
    y2.append(float(split[1][:-1]))
    y3.append(float(split[1][:-1]))
    y4.append(float(split[1][:-1]))


plt.subplot(221)
plt.scatter(x1, y1, s=1)
plt.xlabel('Predicted')
plt.ylabel('Experimental')
plt.title('Solubility from Octanol/Water Partition Coefficient(LogP)')

plt.subplot(222)
plt.scatter(x2, y2, s=1)
plt.xlabel('Predicted')
plt.ylabel('Experimental')
plt.title('Solubility from Atom Pair Fingerprint')

plt.subplot(223)
plt.scatter(x3, y3, s=1)
plt.xlabel('Predicted')
plt.ylabel('Experimental')
plt.title('Solubility from ESOL')

plt.subplot(224)
plt.scatter(x4, y4, s=1)
plt.xlabel('Predicted')
plt.ylabel('Experimental')
plt.title('Solubility From Combined Neural Network')

plt.tight_layout()

plt.show()


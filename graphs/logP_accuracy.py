# Graphs the experimental vs predicted logP
import pickle
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

data = open('data/OctanolWaterPartitionCoefficient/owpc_fixed.txt', 'r')

logP_model = pickle.load(open('run_models/logP_model.pkl', 'rb'))
logP_scaler = pickle.load(open('run_models/logP_scaler.pkl', 'rb'))

x = []
y = []

for line in data.readlines():
    split = line.split(' ')

    compound = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0]))
    logP = logP_model.predict(logP_scaler.transform(np.asarray(compound).reshape(1, -1)))[0]

    x.append(logP)
    y.append(float(split[1][:-1]))

plt.scatter(x, y, s=1)
plt.show()
# Graphs the experimental vs predicted logP
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from common.chemical_models import LogP

data = open('data/OctanolWaterPartitionCoefficient/owpc_fixed.txt', 'r')

# Load necessary model
logP_model = LogP('logP')

x = []
y = []

for line in data.readlines():
    split = line.split(' ')

    # Calculate Atom Pair Fingerprint
    mol = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)

    # Predict logP
    logP = logP_model.run(fingerprint)

    x.append(logP)
    y.append(float(split[1][:-1]))

plt.scatter(x, y, s=1)
plt.show()
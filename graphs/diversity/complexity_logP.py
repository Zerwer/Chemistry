# Produces a 3D scatter plot of the distribution of a list of molecules
#   based on complexity, logP values, and molecular weight
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, GraphDescriptors, Descriptors
from chemical_models import LogP


logP_model = LogP('logP')

smiles = [smile[:-1] for smile in open('data/smiles/all.txt', 'r').readlines()]

x = []
y = []
z = []

new = []

for smile in smiles:
    mol = Chem.MolFromSmiles(smile)
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)

    logP = logP_model.run(fingerprint)
    bertz = GraphDescriptors.BertzCT(mol)

    if logP < 10 and bertz < 1000:
        new.append(smile)
        x.append(logP)
        y.append(bertz)
        z.append(Descriptors.ExactMolWt(mol))


new_file = open('data/smiles/filtered.txt', 'w')
for smile in new:
    new_file.write(smile + '\n')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel('logP')
ax.set_ylabel('Bertz CT')
ax.set_zlabel('Mol Wt')

ax.scatter(x, y, z, s=1, alpha=0.1)
plt.show()


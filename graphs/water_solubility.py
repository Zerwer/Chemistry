# Compares several models ability to predict solubility
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from chemical_models import AtomPairSolubility, LogP, LogPSolubility, CombinedSolubility

data = open('data/water_solubility/aqsol.txt', 'r')

logP_model = LogP('logP')
logP_solubility_model = LogPSolubility('logS_logP')
atom_pair_sol_model = AtomPairSolubility('water_solubility')
combined_model = CombinedSolubility('combined_solubility')

x1, x2, x3, x4 = [], [], [], []
y1, y2, y3, y4 = [], [], [], []

# Read the molecule and corresponding solubility and split into X and Y
for line in data.readlines():
    split = line.split(' ')

    compound = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound)

    logP = logP_model.run(fingerprint)
    logP_sol = logP_solubility_model.run(logP)
    atom_pair_sol = atom_pair_sol_model.run(fingerprint)
    combined_sol = combined_model.run(compound, logP, logP_sol, atom_pair_sol)

    # Additional ESOL empirical model to increase accuracy
    mw = Descriptors.ExactMolWt(compound)
    rb = rdMolDescriptors.CalcNumRotatableBonds(compound)
    ap = len(compound.GetSubstructMatches(Chem.MolFromSmarts('[a]')))/compound.GetNumHeavyAtoms()
    esol = 0.16 - 0.63*logP - 0.0062*mw + 0.066*rb - 0.74*ap

    x1.append(float(logP_sol))
    x2.append(float(atom_pair_sol))
    x3.append(float(esol))
    x4.append(float(combined_sol))

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


"""
Graphs the predicted versus actual melting point using the reverse gse model
"""
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from common.chemical_models import AtomPairSolubility, LogP, LogPSolubility, CombinedSolubility, MeltingPoint

data = open('data/melting_point/mp.txt', 'r')

# Load necessary models
logP_model = LogP('logP')
logP_solubility_model = LogPSolubility('logS_logP')
atom_pair_sol_model = AtomPairSolubility('water_solubility')
combined_model = CombinedSolubility('combined_solubility')
melting_model = MeltingPoint('melting_gse')

x = []
y = []

# Read the molecule and corresponding melting point and split into X and Y
for line in data.readlines():
    split = line.split(' ')

    # Generate RDKit molecule and Atom Pair fingerprint
    compound = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound)

    # Use models to predict logP and logS
    logP = logP_model.run(fingerprint)
    logP_sol = logP_solubility_model.run(logP)
    atom_pair_sol = atom_pair_sol_model.run(fingerprint)
    combined_sol = combined_model.run(compound, logP, logP_sol, atom_pair_sol)

    mp = melting_model.run(combined_sol, logP)

    x.append(float(mp))
    y.append(float(split[1][:-1]))

plt.scatter(x, y, s=1)
plt.xlabel('Predicted')
plt.ylabel('Experimental')
plt.title('Melting Point Prediction')
plt.show()


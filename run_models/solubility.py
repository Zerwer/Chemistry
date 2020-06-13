# Predicts solubility using the combined solubility model
# Usage: python3 solubility.py SMILES
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from chemical_models import AtomPairSolubility, LogP, LogPSolubility, CombinedSolubility

# Load necessary models
logP_model = LogP('logP')
logP_solubility_model = LogPSolubility('logS_logP')
atom_pair_sol_model = AtomPairSolubility('water_solubility')
combined_model = CombinedSolubility('combined_solubility')

# Generate RDKit molecule and Atom Pair fingerprint
compound = Chem.MolFromSmiles(sys.argv[1])
fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound)

# Use models to predict properties
logP = logP_model.run(fingerprint)
logP_sol = logP_solubility_model.run(logP)
atom_pair_sol = atom_pair_sol_model.run(fingerprint)
combined_sol = combined_model.run(compound, logP, logP_sol, atom_pair_sol)

# ESOL estimation from previous predicted data
mw = Descriptors.ExactMolWt(compound)
rb = rdMolDescriptors.CalcNumRotatableBonds(compound)
ap = len(compound.GetSubstructMatches(Chem.MolFromSmarts('[a]'))) / compound.GetNumHeavyAtoms()
esol = 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * rb - 0.74 * ap


# Converts solubility in logS to g/L or mg/mL
def logs_conversion(logS):
    return (10.0**float(logS)) * float(Descriptors.ExactMolWt(Chem.MolFromSmiles(sys.argv[1])))


print('LogP solubility prediction: ' + str(logs_conversion(logP_sol)))
print('Atomic Pair solubility prediction: ' + str(logs_conversion(atom_pair_sol)))
print('ESOL solubility prediction: ' + str(logs_conversion(esol)))
print('Combined solubility prediction: ' + str(logs_conversion(combined_sol)))
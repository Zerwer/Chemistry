# Predicts solubility using the combined solubility model
# Usage: python3 solubility.py SMILES
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from chemical_models import AtomPairSolubility, LogP, LogPSolubility, CombinedSolubility
from functions import logs_to_mg_ml

# Load necessary models
logP_model = LogP('logP')
logP_solubility_model = LogPSolubility('logS_logP')
atom_pair_sol_model = AtomPairSolubility('water_solubility')
combined_model = CombinedSolubility('combined_solubility')

# Generate RDKit molecule and Atom Pair fingerprint
mol = Chem.MolFromSmiles(sys.argv[1])
fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)

# Use models to predict properties
logP = logP_model.run(fingerprint)
logP_sol = logP_solubility_model.run(logP)
atom_pair_sol = atom_pair_sol_model.run(fingerprint)
combined_sol = combined_model.run(mol, logP, logP_sol, atom_pair_sol)

# ESOL estimation from previous predicted data
mw = Descriptors.ExactMolWt(mol)
rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
ap = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[a]'))) / mol.GetNumHeavyAtoms()
esol = 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * rb - 0.74 * ap

# All outputs are converted from logS to mg/mL
print('LogP solubility prediction: ' +
      str(logs_to_mg_ml(logP_sol, mol)) + ' mg/mL')
print('Atomic Pair solubility prediction: ' +
      str(logs_to_mg_ml(atom_pair_sol, mol)) + ' mg/mL')
print('ESOL solubility prediction: ' +
      str(logs_to_mg_ml(esol, mol)) + ' mg/mL')
print('Combined solubility prediction: ' +
      str(logs_to_mg_ml(combined_sol, mol)) + ' mg/mL')

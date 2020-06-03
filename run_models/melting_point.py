# Predicts melting point from logP and solubility
# Usage: python3 melting_point.py SMILES
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from common.chemical_models import AtomPairSolubility, LogP, LogPSolubility, CombinedSolubility, MeltingPoint

# Load necessary models
logP_model = LogP('logP')
logP_solubility_model = LogPSolubility('logS_logP')
atom_pair_sol_model = AtomPairSolubility('water_solubility')
combined_model = CombinedSolubility('combined_solubility')
melting_point_model = MeltingPoint('melting_gse')

# Generate RDKit molecule and Atom Pair fingerprint
compound = Chem.MolFromSmiles(sys.argv[1])
fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound)

# Use models to predict properties
logP = logP_model.run(fingerprint)
logP_sol = logP_solubility_model.run(logP)
atom_pair_sol = atom_pair_sol_model.run(fingerprint)
combined_sol = combined_model.run(compound, logP, logP_sol, atom_pair_sol)
mp = melting_point_model.run(combined_sol, logP)

print(mp)

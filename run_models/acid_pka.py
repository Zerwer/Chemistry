"""
Predict the pKa of an acid from SMILES string
Returns both the fingerprint model prediction and similarity model predictions
"""
from chemical_models import AcidSimilarity, AcidpKa
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, MACCSkeys
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
import sys

# Load models
sim_model = AcidSimilarity('acid_sim')
fp_model = AcidpKa('pKa_acid')

# Set of acids required for similarity model
acid_data = open('data/pKa/formatted_acidic.txt', 'r')
acids = []

mol = Chem.MolFromSmiles(sys.argv[1])

# Read acids from file
for line in acid_data.readlines():
    split = line.split(' ')
    acids.append([split[0], float(split[1][:-1]),
                  rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0]))])

# Run the models and print results
print("Similarity based model: " + str(sim_model.run(sys.argv[1], acids)))
print("Fingerprint based model: " + str(fp_model.run(GetAvalonFP(mol) + MACCSkeys.GenMACCSKeys(mol) +
                                                     rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol))))

# Predicts pKa of any molecule
from chemical_models import GeneralPKa
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, MACCSkeys
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
import sys

# Load models
model = GeneralPKa('pKa')

mol = Chem.MolFromSmiles(sys.argv[1])

# Run the models and print results
print("Predicted pKa: " + str(model.run(GetAvalonFP(mol) +
                                        MACCSkeys.GenMACCSKeys(mol) +
                                        rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol))))
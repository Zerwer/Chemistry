# Usage: python3 run_model.py model_name.pkl scaler_name.pkl SMILE
import pickle, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

model = pickle.load(open(sys.argv[1], 'rb'))
scaler = pickle.load(open(sys.argv[2], 'rb'))

compound = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(sys.argv[3]))
print(model.predict(scaler.transform(np.asarray(compound).reshape(1, -1))))


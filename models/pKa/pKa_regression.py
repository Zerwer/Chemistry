# Predicts pKa by combining three fingerprints, avalon, maacs,
#   and atom pair
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import MACCSkeys
from rdkit.Avalon.pyAvalonTools import GetAvalonFP

acid_data = open('data/pKa/formatted_acidic.txt', 'r')
base_data = open('data/pKa/formatted_basic.txt', 'r')

found_smiles = []
X = []
y = []

# Read the molecule and corresponding pKa and split into X and Y
for line in acid_data.readlines():
    split = line.split(' ')

    mol = Chem.MolFromSmiles(split[0])

    found_smiles.append(split[0])
    # Combine avalon, maacs, atom pair fingerprints
    X.append(GetAvalonFP(mol) +
             MACCSkeys.GenMACCSKeys(mol) +
             rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol))
    y.append(float(split[1][:-1]))

for line in base_data.readlines():
    split = line.split(' ')

    if split[0] in found_smiles:
        continue

    mol = Chem.MolFromSmiles(split[0])

    found_smiles.append(split[0])
    # Combine avalon, maacs, atom pair fingerprints
    X.append(GetAvalonFP(mol) +
             MACCSkeys.GenMACCSKeys(mol) +
             rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol))
    y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
y = np.asarray(y)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1,
                                                    random_state=2)

model = MLPRegressor(solver='adam',
                     alpha=0.0001,
                     hidden_layer_sizes=(512, 256,),
                     random_state=1,
                     verbose=1,
                     max_iter=116,
                     batch_size=500)
model.fit(X_train, y_train)

print(model.score(X_test, y_test))

save_model = open('run_models/pKa_model.pkl', 'wb')
save_scaler = open('run_models/pKa_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)
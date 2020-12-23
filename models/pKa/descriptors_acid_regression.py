# Predict pKa of acids from descriptors and combination of fingerprint
#   and similarity predictors
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Lipinski, MACCSkeys
from chemical_models import AcidpKa, AcidSimilarity
from rdkit.Avalon.pyAvalonTools import GetAvalonFP

data = open('data/pKa/formatted_acidic.txt', 'r')

acids = []

for line in data.readlines():
    split = line.split(' ')
    mol = Chem.MolFromSmiles(split[0])
    acids.append([split[0], float(split[1][:-1]),
                  rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)])

acid_model = AcidpKa('pKa_acid')
sim_model = AcidSimilarity('acid_sim')

X = []
Y = []

# For x combine predictions and descriptors, for y append actual pKa
for line in data.readlines():
    split = line.split(' ')

    mol = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)

    pKa = acid_model.run(GetAvalonFP(mol) +
                         MACCSkeys.GenMACCSKeys(mol) +
                         fingerprint)
    sim_pKa = sim_model.run(split[0], acids)

    X.append([pKa,
              sim_pKa,
              Lipinski.NumHDonors(mol),
              Lipinski.NumHAcceptors(mol),
              Lipinski.NHOHCount(mol)])

    Y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1,
                                                    random_state=1)

model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1048, 128),
                     random_state=1, verbose=1, max_iter=1000, batch_size=500)
model.fit(X_train, y_train)

print(model.score(X_test, y_test))

# Regression model to predict solubility in LogS from the aqsoldb
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

data = open('data/water_solubility/aqsol.txt', 'r')

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')

    mol = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)

    X.append(fingerprint)
    Y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1,
                                                    random_state=1)

model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1026, 128,),
                     random_state=1, verbose=1, max_iter=50, batch_size=100)
model.fit(X_train, y_train)

print(model.score(X_test, y_test))

save_model = open('run_models/water_solubility_model.pkl', 'wb')
save_scaler = open('run_models/water_solubility_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)

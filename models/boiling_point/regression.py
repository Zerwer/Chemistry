"""
Regression neural network to predict boiling point

Requires a deeper network with more samples to achieve any good results

Atom Pair               |   1026, 128
"""
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors

data = open('data/boiling_point/bt.txt', 'r')

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')
    X.append(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))

    # Convert to kelvin
    Y.append(float(split[1][:-1])+273)

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1026, 128,), random_state=1, verbose=1, max_iter=57, batch_size=100)

model.fit(X_train, y_train)

print(model.score(X_test, y_test))
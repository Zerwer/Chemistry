"""
Regression neural network for predicting octanol-water partition coefficient
Used five fingerprints to determine which is best
Result: Atom Pair
layers = (1026, 138), iter = 20, batch = 500
"""
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors

data = open('data/OctanolWaterPartitionCoefficient/owpc_fixed.txt', 'r')

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')
    # X.append(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(split[0]), 2))
    # X.append(Chem.RDKFingerprint(Chem.MolFromSmiles(split[0])))
    # X.append(MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(split[0])))
    X.append(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))
    # X.append(rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))

    Y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

clf = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1026, 128,), random_state=2, verbose=1, max_iter=55, batch_size=500)

clf.fit(X_train, y_train)

print(clf.score(X_test, y_test))

# To predict value for a compound:
# compound = np.asarray(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles('CCO'))).reshape(1, -1)
# print(clf.predict(scaler.transform(compound)))

from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors
data = open('data/OctanolWaterPartitionCoefficient/owpc.txt', 'r')

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

X = preprocessing.scale(np.asarray(X))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

clf = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1026, 128,), random_state=1, verbose=1, max_iter=18, batch_size=10)

clf.fit(X_train, y_train)

print(clf.score(X_test, y_test))
print(clf.predict(np.asarray(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(
    'CC1CCC2CC(C(=CC=CC=CC(CC(C(=O)C(C(C(=CC(C(=O)CC(OC(=O)C3CCCCN3C(=O)C(=O)C1(O2)O)C'
    '(C)CC4CCC(C(C4)OC)OC(=O)C(C)(CO)CO)C)C)O)OC)C)C)C)OC'))).reshape(1, -1)))
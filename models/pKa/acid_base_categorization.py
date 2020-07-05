from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors

acid_data = open('data/pKa/formatted_acidic.txt', 'r')
base_data = open('data/pKa/formatted_basic.txt', 'r')

smiles = []
X = []
y = []

for line in acid_data.readlines():
    split = line.split(' ')
    smiles.append(split[0])
    X.append(MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(split[0])))
    y.append([1, 0, 0])

for line in base_data.readlines():
    split = line.split(' ')

    if split[0] in smiles:
        y.append([0, 0, 1])
    else:
        y.append([0, 1, 0])

    X.append(MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(split[0])))

X = preprocessing.scale(np.asarray(X))
y = np.asarray(y)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=1)

clf = MLPClassifier(solver='adam', alpha=1e-5, hidden_layer_sizes=(96, 32,), random_state=1, verbose=1, max_iter= 2000)

clf.fit(X_train, y_train)

print(clf.score(X_test, y_test))

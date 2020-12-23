# Determines if a molecule will be acidic, basic, or amphoteric using
#   descriptors and MAACS fingerprint
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from rdkit import Chem
from rdkit.Chem import MACCSkeys, Lipinski

# Use combined dataset created using data_collection/baa_sorter.py
data = open('data/pKa/all.txt', 'r')

smiles = []
X = []
y = []

# Read file and set X to relevant molecular data and y to the type, acid, base, or amphoteric
for line in data.readlines():
    split = line.split(' ')
    mol = Chem.MolFromSmiles(split[0])

    X.append(np.append(np.asarray(MACCSkeys.GenMACCSKeys(mol)),
                       [Lipinski.NumHDonors(mol),
                        Lipinski.NumHAcceptors(mol)]))

    if split[1][:-1] == '0':
        y.append([1, 0, 0])
    elif split[1][:-1] == '1':
        y.append([0, 1, 0])
    else:
        y.append([0, 0, 1])

X = preprocessing.scale(np.asarray(X))
y = np.asarray(y)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1,
                                                    random_state=1)

clf = MLPClassifier(solver='adam', alpha=1e-5, hidden_layer_sizes=(128, 32,),
                    random_state=1, verbose=1, max_iter=100)

clf.fit(X_train, y_train)

print(clf.score(X_test, y_test))

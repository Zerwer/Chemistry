# Simple neural network sorter for five fingerprints
# Uncomment selected fingerprints append and set layers based on below
# Purpose was to determine best fingerprint for this application, Atom Pair worked best
#
# Fingerprints and best network shape
# 0 - Morgan                  |   512, 16
# 1 - Daylight                |   513, 64
# 2 - MACCS                   |   96, 32
# 3 - Atom Pair               |   1026, 128
# 4 - Topological Torsion     |   512, 48
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors

fingerprint_type = 3
fingerprints = {
    0: AllChem.GetMorganFingerprintAsBitVect,
    1: Chem.RDKFingerprint,
    2: MACCSkeys.GenMACCSKeys,
    3: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect,
    4: rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect
}
shape = {0: (512, 16),
         1: (513, 64),
         2: (96, 32),
         3: (1026, 128),
         4: (512, 48)}

data = open('data/boiling_point/bt.txt', 'r')

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')
    mol = Chem.MolFromSmiles(split[0])

    kwargs = {0: [mol, 2], 1: [mol], 2: [mol], 3: [mol], 4: [mol]}

    X.append(fingerprints[fingerprint_type](*kwargs[fingerprint_type]))

    # Ranges for boiling point based on normal distribution
    if float(split[1])+273 < 348:
        Y.append([1, 0, 0, 0])
    elif float(split[1])+273 <= 459:
        Y.append([0, 1, 0, 0])
    elif float(split[1])+273 <= 570:
        Y.append([0, 0, 1, 0])
    elif float(split[1])+273 > 570:
        Y.append([0, 0, 0, 1])

X = preprocessing.scale(np.asarray(X))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1,
                                                    random_state=1)

clf = MLPClassifier(solver='adam', alpha=1e-5,
                    hidden_layer_sizes=shape[fingerprint_type],
                    random_state=1, verbose=0, max_iter=2000)
clf.fit(X_train, y_train)

print(clf.score(X_test, y_test))

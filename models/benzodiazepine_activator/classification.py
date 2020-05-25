from sklearn.neural_network import BernoulliRBM
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors

data = open('data/benzodiazepine_activator/total_smiles.txt', 'r')

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')
    # X.append(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(split[0]), 2))
    # X.append(Chem.RDKFingerprint(Chem.MolFromSmiles(split[0])))
    # X.append(MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(split[0])))
    X.append(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))
    # X.append(rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))

    Y.append(1)
    # if float(split[1])+273 < 348:
    #     Y.append([1, 0, 0, 0])
    # elif float(split[1])+273 <= 459:
    #     Y.append([0, 1, 0, 0])
    # elif float(split[1])+273 <= 570:
    #     Y.append([0, 0, 1, 0])
    # elif float(split[1])+273 > 570:
    #     Y.append([0, 0, 0, 1])

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

model = BernoulliRBM(random_state=1, verbose=1)

model.fit(X_train)

print(model.score_samples(X_test))
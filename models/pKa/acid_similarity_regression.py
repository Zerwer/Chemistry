import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors, Lipinski
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn import preprocessing


def five_closest(mol, acid_set):
    mol_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(mol))
    similarity = []
    for acid in acid_set:
        sim = DataStructs.DiceSimilarity(mol_fp, acid[2])
        similarity.append([sim, acid[1]])

    return np.asarray(sorted(similarity)[:5]).flatten()

acid_data = open('data/pKa/formatted_acidic.txt', 'r')
acids = []

for line in acid_data.readlines():
    split = line.split(' ')
    acids.append([split[0], float(split[1][:-1]), rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0]))])

# Split data into training and test set
train, test = train_test_split(acids, test_size=0.25, random_state=1)

X = []
y = []

for acid in test:
    X.append(five_closest(acid[0], acids))
    y.append(acid[1])

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
y = np.asarray(y)

# Split data into training and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

model = MLPRegressor(solver='adam',
                     alpha=0.0001,
                     hidden_layer_sizes=(256, 128,),
                     random_state=1,
                     verbose=1,
                     max_iter=1000)
model.fit(X, y)
print(model.score(X_test, y_test))
# Regression neural network for predicting octanol-water partition coefficient
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

data = open('data/OctanolWaterPartitionCoefficient/owpc_fixed.txt', 'r')

X = []
y = []

# Read the molecule and corresponding logP and split into X and Y
for line in data.readlines():
    split = line.split(' ')

    mol = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
    X.append(fingerprint)
    y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
y = np.asarray(y)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=1)

model = MLPRegressor(solver='adam',
                     alpha=0.0001,
                     hidden_layer_sizes=(1026, 128,),
                     random_state=0,
                     verbose=1,
                     max_iter=200,
                     batch_size=500)
model.fit(X_train, y_train)

print(model.score(X_test, y_test))

save_model = open('run_models/logP_model.pkl', 'wb')
save_scaler = open('run_models/logP_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)

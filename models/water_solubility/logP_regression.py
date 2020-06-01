"""
Regression prediction of LogS from Octanol-Water partition coefficient (LogP)
"""
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from models.common.chemical_models import LogP

data = open('data/water_solubility/aqsol.txt', 'r')

logP_model = LogP('logP')

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')

    compound = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0]))

    logP = logP_model.run(compound)

    X.append(logP)

    Y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X).reshape(-1, 1))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1048, 256), random_state=1, verbose=1, max_iter=34, batch_size=500)

model.fit(X_train, y_train)

print(model.score(X_test, y_test))

save_model = open('models/logS_logP_model.pkl', 'wb')
save_scaler = open('models/logS_logP_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)
"""
This model predicts solubility of a molecule by combining the logP regression model and Fingerprint regression model
"""
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

data = open('data/water_solubility/aqsol.txt', 'r')

logP_model = pickle.load(open('models/logP_model.pkl', 'rb'))
logP_scaler = pickle.load(open('models/logP_scaler.pkl', 'rb'))

logP_solubility_model = pickle.load(open('models/logS_logP_model.pkl', 'rb'))
logP_solubility_scaler = pickle.load(open('models/logS_logP_scaler.pkl', 'rb'))

solubility_model = pickle.load(open('models/water_solubility_model.pkl', 'rb'))
solubility_scaler = pickle.load(open('models/water_solubility_scaler.pkl', 'rb'))

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')
    compound = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0]))
    logP = logP_model.predict(logP_scaler.transform(np.asarray(compound).reshape(1, -1)))[0]
    logP_sol = logP_solubility_model.predict(logP_solubility_scaler.transform(np.asarray(logP).reshape(1, -1)))[0]
    sol = solubility_model.predict(solubility_scaler.transform(np.asarray(compound).reshape(1, -1)))[0]

    X.append([logP_sol, sol])

    Y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1048, 128), random_state=1, verbose=1, max_iter=25, batch_size=500)
model.fit(X_train, y_train)

print(model.score(X_test, y_test))

save_model = open('models/combined_solubility_model.pkl', 'wb')
save_scaler = open('models/combined_solubility_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)
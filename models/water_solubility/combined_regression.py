"""
This model predicts solubility of a molecule by combining the logP regression model and Fingerprint regression model
"""
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors

data = open('data/water_solubility/aqsol.txt', 'r')

logP_model = pickle.load(open('run_models/logP_model.pkl', 'rb'))
logP_scaler = pickle.load(open('run_models/logP_scaler.pkl', 'rb'))

logP_solubility_model = pickle.load(open('run_models/logS_logP_model.pkl', 'rb'))
logP_solubility_scaler = pickle.load(open('run_models/logS_logP_scaler.pkl', 'rb'))

solubility_model = pickle.load(open('run_models/water_solubility_model.pkl', 'rb'))
solubility_scaler = pickle.load(open('run_models/water_solubility_scaler.pkl', 'rb'))

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')
    compound = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound)
    logP = logP_model.predict(logP_scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]
    logP_sol = logP_solubility_model.predict(logP_solubility_scaler.transform(np.asarray(logP).reshape(1, -1)))[0]
    sol = solubility_model.predict(solubility_scaler.transform(np.asarray(fingerprint).reshape(1, -1)))[0]

    # Additional ESOL empirical model to increase accuracy
    mw = Descriptors.ExactMolWt(compound)
    rb = rdMolDescriptors.CalcNumRotatableBonds(compound)
    ap = len(compound.GetSubstructMatches(Chem.MolFromSmarts('[a]')))/compound.GetNumHeavyAtoms()
    esol = 0.16 -  0.63*logP - 0.0062*mw + 0.066*rb - 0.74*ap

    X.append([logP_sol, sol, esol])

    Y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1048, 128), random_state=1, verbose=1, max_iter=42, batch_size=500)
model.fit(X_train, y_train)

print(model.score(X_test, y_test))

save_model = open('run_models/combined_solubility_model.pkl', 'wb')
save_scaler = open('run_models/combined_solubility_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)
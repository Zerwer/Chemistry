# This model predicts solubility of a molecule by combining the logP
#   regression model and Fingerprint regression model
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from chemical_models import AtomPairSolubility, LogP, LogPSolubility

data = open('data/water_solubility/aqsol.txt', 'r')

logP_model = LogP('logP')
logP_solubility_model = LogPSolubility('logS_logP')
atom_pair_sol_model = AtomPairSolubility('water_solubility')

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')

    mol = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)

    logP = logP_model.run(fingerprint)
    logP_sol = logP_solubility_model.run(logP)
    atom_pair_sol = atom_pair_sol_model.run(fingerprint)

    # Additional ESOL empirical model to increase accuracy
    mw = Descriptors.ExactMolWt(mol)
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    ap = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[a]'))) / mol.GetNumHeavyAtoms()
    esol = 0.16 - 0.63*logP - 0.0062*mw + 0.066*rb - 0.74*ap

    X.append([logP_sol, atom_pair_sol, esol])
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
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from models.common import models

data = open('data/boiling_point/ia_data.txt', 'r')

logP_model = models.init('logP')
logP_solubility_model = models.init('logS_logP')
atom_pair_sol_model = models.init('water_solubility')
combined_model = models.init('combined_solubility')

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')

    compound = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound)

    logP = models.run_logP(logP_model, fingerprint)
    logP_sol = models.run_logP_sol(logP_solubility_model, logP)
    atom_pair_sol = models.run_atom_pair_solubility(atom_pair_sol_model, fingerprint)
    combined_sol = models.run_combined_solubility(combined_model, compound, logP, logP_sol, atom_pair_sol)

    X.append([combined_sol, logP])

    Y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1048, 128), random_state=1, verbose=1, max_iter=1000, batch_size=500)
model.fit(X_train, y_train)

print(model.score(X_test, y_test))

save_model = open('run_models/boiling_gse_model.pkl', 'wb')
save_scaler = open('run_models/boiling_gse_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)
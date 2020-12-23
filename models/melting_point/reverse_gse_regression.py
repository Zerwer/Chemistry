# Calculates melting point from predicted solubility and octanol water
#   partition coefficient. This is the reverse process of the general
#   solubility equation as it is harder to predict melting point then
#   logS and logP
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from chemical_models import AtomPairSolubility, LogP, LogPSolubility, CombinedSolubility

data = open('data/melting_point/mp.txt', 'r')

logP_model = LogP('logP')
logP_solubility_model = LogPSolubility('logS_logP')
atom_pair_sol_model = AtomPairSolubility('water_solubility')
combined_model = CombinedSolubility('combined_solubility')

X = []
Y = []

# Read the molecule and corresponding solubility and split into X and Y
# X is predicted solubility(logS) and predicted octanol water partition coefficient(logP)
# Y is the melting point
for line in data.readlines():
    split = line.split(' ')

    mol = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)

    logP = logP_model.run(fingerprint)
    logP_sol = logP_solubility_model.run(logP)
    atom_pair_sol = atom_pair_sol_model.run(fingerprint)
    combined_sol = combined_model.run(mol, logP, logP_sol, atom_pair_sol)

    X.append([combined_sol, logP])
    Y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1048, 128), random_state=1, verbose=1, max_iter=1000, batch_size=500)
model.fit(X_train, y_train)

print(model.score(X_test, y_test))

save_model = open('run_models/melting_gse_model.pkl', 'wb')
save_scaler = open('run_models/melting_gse_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)

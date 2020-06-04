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
from common.chemical_models import AtomPairSolubility, LogP, LogPSolubility

data = open('data/water_solubility/aqsol.txt', 'r')

# Load necessary models
logP_model = LogP('logP')
logP_solubility_model = LogPSolubility('logS_logP')
atom_pair_sol_model = AtomPairSolubility('water_solubility')

X = []
Y = []

# Read the molecule and corresponding solubility and split into X and Y
# X is data relevant to the molecule
# Y is the solubility
for line in data.readlines():
    split = line.split(' ')

    # Generate RDKit molecule and Atom Pair fingerprint
    compound = Chem.MolFromSmiles(split[0])
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound)

    # Use models to predict relevant properties through different methods
    logP = logP_model.run(fingerprint)
    logP_sol = logP_solubility_model.run(logP)
    atom_pair_sol = atom_pair_sol_model.run(fingerprint)

    # Additional ESOL empirical model to increase accuracy
    mw = Descriptors.ExactMolWt(compound)
    rb = rdMolDescriptors.CalcNumRotatableBonds(compound)
    ap = len(compound.GetSubstructMatches(Chem.MolFromSmarts('[a]')))/compound.GetNumHeavyAtoms()
    esol = 0.16 - 0.63*logP - 0.0062*mw + 0.066*rb - 0.74*ap

    X.append([logP_sol, atom_pair_sol, esol])
    Y.append(float(split[1][:-1]))

# Scale data for better accuracy
scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

# Split data into training and test set
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

# Create the model and fit to data
model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1048, 128), random_state=1, verbose=1, max_iter=42, batch_size=500)
model.fit(X_train, y_train)

# Score the model off test data that it was not trained on
print(model.score(X_test, y_test))

# Save model and scaler
save_model = open('run_models/combined_solubility_model.pkl', 'wb')
save_scaler = open('run_models/combined_solubility_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)
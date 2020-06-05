"""
Regression model to predict solubility in LogS from the aqsoldb
Atom Pair, (1026, 128), itr = 50, bt = 100, r^2=0.7637
"""
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

data = open('data/water_solubility/aqsol.txt', 'r')

X = []
Y = []

# Read the molecule and corresponding solubility and split into X and Y
for line in data.readlines():
    split = line.split(' ')

    # Use Atom Pair fingerprint as input
    X.append(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))
    Y.append(float(split[1][:-1]))

# Scale data for better accuracy
scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

# Split data into training and test set
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

# Create the model and fit to data
model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1026, 128,), random_state=1, verbose=1, max_iter=50, batch_size=100)
model.fit(X_train, y_train)

# Score the model off test data that it was not trained on
print(model.score(X_test, y_test))

# Save model and scaler
save_model = open('run_models/water_solubility_model.pkl', 'wb')
save_scaler = open('run_models/water_solubility_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)

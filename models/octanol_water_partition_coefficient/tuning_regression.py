from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import GridSearchCV
from sklearn import preprocessing
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

data = open('data/OctanolWaterPartitionCoefficient/owpc_fixed.txt', 'r')

X = []
y = []

# Read the molecule and corresponding logP and split into X and Y
for line in data.readlines():
    split = line.split(' ')

    # Use Atom Pair fingerprint to represent the molecule
    X.append(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))
    y.append(float(split[1][:-1]))

# Scale data for better accuracy
scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
y = np.asarray(y)

# Grid search for optimal hyper-parameters
param_grid = {'alpha': 10.0 ** -np.arange(1, 7),
              'max_iter': [50, 100, 150],
              'batch_size': [500, 3000, 6000]}

# Create the model
model = MLPRegressor(solver='adam', hidden_layer_sizes=(1026, 128,), verbose=1)

# Using Grid Search and Cross Validation find best parameters, use all available cores
search = GridSearchCV(model, param_grid=param_grid, cv=5, n_jobs=-1)
search.fit(X, y)

# Print the best parameters
print(search.best_params_)

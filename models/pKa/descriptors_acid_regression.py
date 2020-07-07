"""
Predict pKa of acids from descriptors and combination of fingerprint and similarity predictors
"""
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Lipinski, MACCSkeys
from chemical_models import AcidpKa, AcidSimilarity
from rdkit.Avalon.pyAvalonTools import GetAvalonFP

data = open('data/pKa/formatted_acidic.txt', 'r')

acids = []

# Read molecules and add to acids for prediction reference
for line in data.readlines():
    split = line.split(' ')
    acids.append([split[0], float(split[1][:-1]),
                  rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0]))])

# Load necessary models
acid_model = AcidpKa('pKa_acid')
sim_model = AcidSimilarity('acid_sim')

X = []
Y = []

# For x combine predictions and descriptors, for y append actual pKa
for line in data.readlines():
    split = line.split(' ')

    compound = Chem.MolFromSmiles(split[0])

    pKa = acid_model.run(GetAvalonFP(compound)+MACCSkeys.GenMACCSKeys(compound) +
                         rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound))
    sim_pKa = sim_model.run(split[0], acids)

    X.append([pKa, sim_pKa, Lipinski.NumHDonors(compound), Lipinski.NumHAcceptors(compound), Lipinski.NHOHCount(compound)])

    Y.append(float(split[1][:-1]))

# Scale data for better accuracy
scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

# Split data into training and test set
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=1)

# Create the model and fit to data
model = MLPRegressor(solver='adam', alpha=1e-5, hidden_layer_sizes=(1048, 128), random_state=1, verbose=1, max_iter=1000, batch_size=500)
model.fit(X_train, y_train)

# Score the model off test data that it was not trained on
print(model.score(X_test, y_test))

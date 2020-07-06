from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Lipinski, MACCSkeys
from chemical_models import AcidpKa
from rdkit.Avalon.pyAvalonTools import GetAvalonFP

data = open('data/pKa/formatted_acidic.txt', 'r')

# Load necessary models
acid_model = AcidpKa('pKa_acid')

X = []
Y = []

for line in data.readlines():
    split = line.split(' ')

    compound = Chem.MolFromSmiles(split[0])

    pKa = acid_model.run(GetAvalonFP(compound)+MACCSkeys.GenMACCSKeys(compound) +
                         rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound))

    X.append([pKa, Lipinski.NumHDonors(compound), Lipinski.NumHAcceptors(compound), Lipinski.NHOHCount(compound)])

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

# Save model and scaler
# save_model = open('run_models/boiling_gse_model.pkl', 'wb')
# save_scaler = open('run_models/boiling_gse_scaler.pkl', 'wb')
# pickle.dump(model, save_model)
# pickle.dump(scaler, save_scaler)
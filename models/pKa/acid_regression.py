from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn import preprocessing
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Lipinski
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit.Avalon.pyAvalonTools import GetAvalonFP

data = open('data/pKa/formatted_acidic.txt', 'r')

X = []
y = []

# Read the molecule and corresponding pKa and split into X and Y
for line in data.readlines():
    split = line.split(' ')

    # Use Atom Pair fingerprint to represent the molecule
    X.append(GetAvalonFP(Chem.MolFromSmiles(split[0]))+MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(split[0])) +
             rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))
    #X.append(MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(split[0])))
    #X.append(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))
    y.append(float(split[1][:-1]))

# Scale data for better accuracy
scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
y = np.asarray(y)

# Split data into training and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=1)

# Create the model and fit to data
model = MLPRegressor(solver='adam',
                     alpha=0.0001,
                     hidden_layer_sizes=(256, 128,),
                     random_state=1,
                     verbose=1,
                     max_iter=116,
                     batch_size=500)
model.fit(X_train, y_train)

# Score the model off test data that it was not trained on
print(model.score(X_test, y_test))

# Save model and scaler
save_model = open('run_models/pKa_acid_model.pkl', 'wb')
save_scaler = open('run_models/pKa_acid_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)
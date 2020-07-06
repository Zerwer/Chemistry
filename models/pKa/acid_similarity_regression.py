import numpy as np
import pickle
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors, Lipinski
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn import preprocessing
from chemical_models import pka_similarities

acid_data = open('data/pKa/formatted_acidic.txt', 'r')
acids = []

for line in acid_data.readlines():
    split = line.split(' ')
    acids.append([split[0], float(split[1][:-1]),
                  rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0]))])

# Split data into training and test set
train, test = train_test_split(acids, test_size=0.5, random_state=1)

X = []
y = []

for acid in test:
    X.append(pka_similarities(acid[0], train, 512))
    y.append(acid[1])

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
y = np.asarray(y)

# Split data into training and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

model = MLPRegressor(solver='adam',
                     alpha=0.0001,
                     hidden_layer_sizes=(512, 128,),
                     random_state=1,
                     verbose=1,
                     max_iter=1000)
model.fit(X_train, y_train)
print(model.score(X_test, y_test))

# Save model and scaler
save_model = open('run_models/acid_sim_model.pkl', 'wb')
save_scaler = open('run_models/acid_sim_scaler.pkl', 'wb')
pickle.dump(model, save_model)
pickle.dump(scaler, save_scaler)
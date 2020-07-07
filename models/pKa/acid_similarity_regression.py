"""
pKa regression that uses similar molecules with known pKa values as input
Either works very well for some molecules if similar ones are known and very poorly for others when they are not
Due to lack of data and the majority of the split going to reference the accuracy is artificially low
"""
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn import preprocessing
from functions import pka_similarities

acid_data = open('data/pKa/formatted_acidic.txt', 'r')
acids = []

# Read acids from data and store SMILES, pKa, fingerprint
for line in acid_data.readlines():
    split = line.split(' ')
    # Atom pair fingerprint used to determine similarity
    acids.append([split[0], float(split[1][:-1]),
                  rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0]))])

# Split data into reference set that will be used to get similarity and
# test set which will be used to train and validate the model
reference, test = train_test_split(acids, test_size=0.5, random_state=1)

X = []
y = []

# X is 512 of the most similar molecules see common.functions.pka_similarities
for acid in test:
    X.append(pka_similarities(acid[0], reference, 512))
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
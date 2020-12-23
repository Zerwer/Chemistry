# Deep neural network for prediction logP values using TensorFlow and Keras
from keras.layers import Dense, Dropout
from keras.models import Sequential
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sklearn import preprocessing

data = open('data/OctanolWaterPartitionCoefficient/owpc_fixed.txt', 'r')

X = []
Y = []

# Read the molecule and corresponding logP and split into X and Y
for line in data.readlines():
    split = line.split(' ')

    X.append(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))
    Y.append(float(split[1][:-1]))

scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

model = Sequential([
    Dense(1026, input_dim=2048,activation='relu'),
    Dropout(0.5),
    Dense(512, activation='relu'),
    Dense(128, activation='relu'),
    Dropout(0.5),
    Dense(1, activation='linear')
])

model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mean_squared_error'])

model.fit(X, Y, validation_split=0.1, shuffle=True, batch_size=500, epochs=100)

print(model.evaluate(X, Y))
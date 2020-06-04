from keras.layers import Dense, Flatten, Dropout
from keras.models import Sequential
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sklearn import preprocessing
from keras import backend as K

data = open('data/OctanolWaterPartitionCoefficient/owpc_fixed.txt', 'r')

X = []
Y = []


def coeff_determination(y_true, y_pred):
    SS_res = K.sum(K.square(y_true-y_pred))
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )

# Read the molecule and corresponding logP and split into X and Y
for line in data.readlines():
    split = line.split(' ')

    # Use Atom Pair fingerprint to represent the molecule
    X.append(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0])))
    Y.append(float(split[1][:-1]))

# Scale data for better accuracy
scaler = preprocessing.StandardScaler()
X = scaler.fit_transform(np.asarray(X))
Y = np.asarray(Y)

model = Sequential([
    Dense(1026, input_dim=2048,activation='relu'),
    Dropout(0.5),
    Dense(128, activation='relu'),
    Dropout(0.5),
    Dense(1, activation='linear')
])

model.compile(optimizer='adam', loss='mean_squared_error', metrics=[coeff_determination])

model.fit(X, Y, validation_split=0.1, shuffle=True, batch_size=500, epochs=100)
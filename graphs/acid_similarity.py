"""
Graph to visualize predicted versus actual acid pKa for similarity model
Splits the same way as the model training
"""
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from chemical_models import AcidSimilarity
from sklearn.model_selection import train_test_split

# Load model
acid_model = AcidSimilarity('acid_sim')

acid_data = open('data/pKa/formatted_acidic.txt', 'r')
acids = []

# Read file to gather reference acids
for line in acid_data.readlines():
    split = line.split(' ')
    acids.append([split[0], float(split[1][:-1]),
                  rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0]))])

# Split data into reference set that will be used to get similarity and
# test set which will be used to train and validate the model
reference, test = train_test_split(acids, test_size=0.5, random_state=1)

X = []
y = []

# Set x to predicted values and y to actual
for acid in test:
    X.append(acid_model.run(acid[0], acids))
    y.append(acid[1])

# Split data into training and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

# Plot the data used to train the model as red and validation data as blue
# If overfit the red data will fit the line significantly more frequently than blue
plt.scatter(X_train, y_train, s=1, color='red')
plt.scatter(X_test, y_test, s=1, color='blue')
plt.show()

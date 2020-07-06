import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from chemical_models import AcidSimilarity
from sklearn.model_selection import train_test_split

acid_model = AcidSimilarity('acid_sim')

acid_data = open('data/pKa/formatted_acidic.txt', 'r')
acids = []

for line in acid_data.readlines():
    split = line.split(' ')
    acids.append([split[0], float(split[1][:-1]), rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(split[0]))])

# Split data into training and test set
reference, test = train_test_split(acids, test_size=0.25, random_state=1)

X = []
y = []

for acid in test:
    X.append(acid_model.run(acid[0], acids))
    y.append(acid[1])

# Split data into training and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

plt.scatter(X_train, y_train, s=1, color='red')
plt.scatter(X_test, y_test, s=1, color='blue')
plt.show()

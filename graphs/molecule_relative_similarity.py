# Creates a png of a list of SMILES
# Iterates over compounds, finds most similar then finds most similar to that
#   compound until all compounds are arranged
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors


def mol_similarity_grid(mols, size, row):
    sort = sort_similarity(mols[1:], [mols[0]])
    sorted_mols = [x[1] for x in sort]
    draw = Draw.MolsToGridImage(sorted_mols, molsPerRow=row, subImgSize=size)
    return draw


# Sorts a list of molecules by similarity using simple insertion sort
def sort_similarity(mols, sort):
    largest = [0, mols[0]]

    # Find the most similar molecule by finding largest similarity score
    for mol in mols:
        similarity = DataStructs.DiceSimilarity(sort[len(sort)-1][0], mol[0])
        if similarity > largest[0]:
            largest = [similarity, mol]

    # Move the molecule from unsorted list to sorted list
    mols.remove(largest[1])
    sort.append(largest[1])

    # check if there are more mols to sort
    if len(mols) > 0:
        sort_similarity(mols, sort)
    return sort


if __name__ == "__main__":
    data = open('data/benzodiazepine_activator/total_smiles.txt', 'r')
    molecules = []

    for line in data.readlines():
        compound = Chem.MolFromSmiles(line[:-1])
        molecules.append((rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(compound), compound))

    new_mols = sort_similarity(molecules[1:], [molecules[0]])
    img = Draw.MolsToGridImage([x[1] for x in new_mols], molsPerRow=10, subImgSize=(200, 200))
    img.save('molecules.png')

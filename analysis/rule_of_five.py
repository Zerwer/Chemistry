# Function that determines if a molecule satisfies Lipinski's rule of five:
# * No more than five hydrogen bond donors (NH and OH bonds)
# * No more than ten hydrogen bond acceptors (N and O atoms)
# * Molecular mass less than 500 daltons
# * LogP value less than 5
import sys
from chemical_models import LogP
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski


def lipinski_rule(mol):
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
    return [
        Lipinski.NHOHCount(mol) <= 5,
        Lipinski.NOCount(mol) <= 10,
        Descriptors.ExactMolWt(mol) <= 500,
        LogP('logP').run(fingerprint) <= 5]


conditions = lipinski_rule(Chem.MolFromSmiles(sys.argv[1]))
print(
    'Less than five hydrogen bond donors: ' + str(conditions[0]) +
    '\nLess than ten hydrogen bond acceptors: ' + str(conditions[1]) +
    '\nMolecular mass smaller than 500: ' + str(conditions[2]) +
    '\nOctanol-water partition coefficient less than 5: ' + str(conditions[3]))

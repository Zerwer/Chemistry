import sys
from chemical_models import LogP
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski


def lipinski_rule(mol, logP_model):
    return [
        Lipinski.NHOHCount(mol) <= 5,
        Lipinski.NOCount(mol) <= 10,
        Descriptors.ExactMolWt(mol) <= 500,
        logP_model.run(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)) <= 5]


conditions = lipinski_rule(Chem.MolFromSmiles(sys.argv[1]), LogP('logP'))
print('Less than five hydrogen bond donors: ' + str(conditions[0]) +
      '\nLess than ten hydrogen bond acceptors: ' + str(conditions[1]) +
      '\nMolecular mass smaller than 500: ' + str(conditions[2]) +
      '\nOctanol-water partition coefficient less than 5: ' + str(conditions[3]))

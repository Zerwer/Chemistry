# Model that predicts structure from fingerprint


def smiles_to_num(smiles):
    return [ord(c) for c in smiles]


def num_to_smiles(num):
    return ''.join([chr(i) for i in num])


print(smiles_to_num("CCO"))
print(num_to_smiles(smiles_to_num("CCO")))
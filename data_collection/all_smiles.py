files = ['data/benzodiazepine_activator/total_smiles.txt',
         'data/boiling_point/all_compounds.txt',
         'data/melting_point/mp.txt',
         'data/OctanolWaterPartitionCoefficient/owpc_fixed.txt',
         'data/pKa/all.txt',
         'data/water_solubility/bkws.txt']

smiles = []

for file in files:
    for line in open(file, 'r').readlines():
        smile = line.split(' ')[0]
        if smile not in smiles:
            smiles.append(smile)


new_file = open('data/smiles/all.txt', 'w')
for smile in smiles:
    new_file.write(smile + '\n')

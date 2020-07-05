file = open('data/pKa/pKaInWater.dwar', 'r')
new_acid = open('data/pKa/formated_acidic.txt', 'w')
new_basic = open('data/pKa/formated_basic.txt', 'w')

for line in file.readlines():
    columns = line.split('\t')
    smiles = columns[len(columns) - 1][:-1]
    pka = columns[4]

    if columns[9] == 'acidic':
        new_acid.write(smiles + ' ' + pka + '\n')
    else:
        new_basic.write(smiles + ' ' + pka + '\n')

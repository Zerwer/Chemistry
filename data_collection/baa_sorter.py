# Combines both acidic and basic data into one file that is either basic,
#   acidic or amphoteric
acid_data = open('data/pKa/formatted_acidic.txt', 'r')
base_data = open('data/pKa/formatted_basic.txt', 'r')

new_data = open('data/pKa/all.txt', 'w')

acids = []
bases = []
amphoteric = []

# Gather all acids and prevent duplicates
for line in acid_data.readlines():
    split = line.split(' ')
    if split[0] not in acids:
        acids.append(split[0])

# Gather all bases, if already in acids remove and add to amphoteric,
#   prevent duplicates for both bases and amphoterics
for line in base_data.readlines():
    split = line.split(' ')
    if split[0] in acids:
        acids.remove(split[0])
        amphoteric.append(split[0])
    elif split[0] not in amphoteric and split[0] not in bases:
        bases.append(split[0])

for acid in acids:
    new_data.write(acid + ' 0\n')
for base in bases:
    new_data.write(base + ' 1\n')
for both in amphoteric:
    new_data.write(both + ' 2\n')

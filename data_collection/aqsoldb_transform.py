aqsoldb = open('data/water_solubility/aqsol.csv', 'r')
new = open('aqsol.txt', 'w')

aqsoldb.readline()

for line in aqsoldb.readlines():
    split = line.split(',')
    new.write(split[len(split) - 2] + ' ' + split[len(split) - 1])
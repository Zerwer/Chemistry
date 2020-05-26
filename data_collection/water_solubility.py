from pubchem import experimental_properties
from collect import collect
import re


def water_solubility(data):
    part_data = experimental_properties(data, 'Solubility')

    # Must be in mg/mL or g/L
    solubilities = []

    for point in part_data['Information']:
        line = str(point['Value']['StringWithMarkup'][0]['String']).lower()
        # Scientific notation g/L
        if re.match(r"^-?\d*\.?\d+e[+-]?\d+ g/l$", line):
            solubilities.append(float(line[:-4]))
        # Decimal g/L
        elif re.match(r"^-?\d*\.?\d+ g/l$", line):
            solubilities.append(float(line[:-4]))
        # Scientific notation mg/mL
        elif re.match(r"^-?\d*\.?\d+e[+-]?\d+ mg/ml$", line):
            solubilities.append(float(line[:-6]))
        # Decimal mg/mL
        elif re.match(r"^-?\d*\.?\d+ mg/ml$", line):
            solubilities.append(float(line[:-6]))
        # Scientific notation mg/L
        elif re.match(r"^-?\d*\.?\d+e[+-]?\d+ mg/l$", line):
            solubilities.append(float(line[:-5])/1000)
        # Decimal mg/L
        elif re.match(r"^-?\d*\.?\d+ mg/l$", line):
            solubilities.append(float(line[:-5])/1000)
        # Decimal mg/mL explicit temperature
        elif re.match(r"^-?\d*\.?\d+ mg/ml at 25 °c$", line):
            solubilities.append(float(line[:-9]))
        # Written scientific mg/L explicit temperature
        elif re.match(r"^in water, -?\d*\.?\d+x10[+-]?\d+ mg/l at 2\d °c$", line):
            solubilities.append(float(line.replace('x10', 'e')[10:-14])/1000)
        # Written scientific mg/L explicit temperature no comma
        elif re.match(r"^in water -?\d*\.?\d+x10[+-]?\d+ mg/l at 2\d °c$", line):
            solubilities.append(float(line.replace('x10', 'e')[9:-14])/1000)
        elif re.match(r"^in water, \d*[,.]?\d+ mg/ml at 2\d °c$", line):
            solubilities.append(float(line[10:-15]))
        # Written float mg/L comma instead of . and @ instead of "at"
        elif re.match(r"^in water, -?\d*,?\d+x10[+-]?\d+ mg/l @ 25 °c$", line):
            solubilities.append(float(line.replace('x10', 'e')[10:-13])/1000)
        # Written integer mg/L
        elif re.match(r"^in water, \d*[,.]?\d+ mg/l at 2\d °c$", line):
            solubilities.append(float(line.replace(',', '.')[10:-14])/1000)
        # Written integer mg/L with ph
        elif re.match(r"^in water, \d*[,.]?\d+ mg/l at 2\d °c, ph 7$", line):
            solubilities.append(float(line.replace(',', '.')[10:-20])/1000)
        # Written long text decimal
        elif re.match(r"^solubility in water, g/100ml at 20 °c: \d*\.\d+$", line):
            solubilities.append(float(line[39:])*10)
        elif re.match(r"^\d*\.?\d+ \[ug\/ml\]$", line):
            solubilities.append(float(line[:-8])/1000)
        elif re.match(r"^[<>]?\d*\.?\d+ \[ug\/ml\]$", line):
            solubilities.append(float(line[1:-8]) / 1000)
        elif re.match(r"^h2o \d*\.?\d+ \( ?mg\/ml\)$", line):
            solubilities.append(float(line.replace(' ', '')[3:-7]))
        # Insoluble
        elif 'water' in line and ('insoluble' in line or 'none' in line) or line == 'insoluble':
            solubilities.append(0)
        # elif len(line) < 20 or ('water' in line and 'est' not in line):
        #     print(line)
        # else:
        #     print(line)

    return round(sum(solubilities)/len(solubilities), 4)

# Takes list of molecules and outputs the smiles and corresponding solubility in water of a compound
collect('data/water_solubility/list.txt', 'water_solubility.txt', water_solubility)
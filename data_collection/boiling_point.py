# Scraping tool to collect boiling point data from PubChem
from pubchem import experimental_properties
from collect import collect, APIFunction
import re


class BoilingPoint(APIFunction):

    # Takes JSON data and finds boilings points and averages them
    def run(self, data):
        boiling_data = str(experimental_properties(data, 'Boiling Point'))

        bp = []

        # Finds each boiling point and uses different regular expression
        #   searches based on format
        for line in boiling_data.split('\''):
            if '°C' in line and line != '°C':
                if re.search(r'>\d+', line) or re.search(r'<\d+', line):
                    pass
                elif '±' in line:
                    bp.append(float(re.findall(r'\d+.?\d+±', line)[0][:-1]))
                elif '@' in line or 'kPa' in line:
                    print(line)
                elif re.search(r'[0-9]-[0-9]', line):
                    x = re.findall(r'\d+.?\d+-\d+.?\d+', line)[0].split('-')
                    bp.append((float(x[0])+float(x[1]))/2)
                else:
                    bp.append(float(re.sub('[ °C]', '',
                            re.findall(r'-?\d+(?:\.\d+)?\s*°\s*C', line)[0])))

        # Prevent division by zero
        if len(bp) == 0:
            return None
        else:
            return sum(bp)/len(bp)


# Takes list of molecules and outputs the smiles and corresponding boiling
#   point from boiling_point function
collect('list.txt', 'bp.txt', BoilingPoint)

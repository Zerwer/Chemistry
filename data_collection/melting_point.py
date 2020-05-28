from pubchem import experimental_properties
from collect import collect
import re


# Takes JSON data and finds melting points and averages them
def melting_point(data):
    melting_data = experimental_properties(data, 'Melting Point')

    mp = []

    for point in melting_data['Information']:
        line = str(point['Value']['StringWithMarkup'][0]['String']).lower()
        if 'decomp' not in line and \
           'greater' not in line and \
           'less' not in line and \
           'approx' not in line and \
           '<' not in line and \
           '>' not in line:

            # Get either a range of temperature or the measured temperature
            if '°c' in line:
                if re.match(r"\d*-\d+", line):
                    ranges = re.findall(r"\d*-\d+", line)[0].split['-']
                    mp.append((float(ranges[0]) + float(ranges[1]))/2)
                else:
                    mp.append(float(re.findall(r"\d*", line)[0]))
            # Same as getting celsius just convert
            elif '°f' in line:
                if re.match(r"\d*-\d+", line):
                    ranges = re.findall(r"\d*-\d+", line)[0].split['-']
                    mp.append((((float(ranges[0]) + float(ranges[1]))/2)-32) * (5/9))
                else:
                    mp.append((float(re.findall(r"\d*", line)[0]) - 32) * (5/9))

    if len(mp) == 0:
        return None
    else:
        return sum(mp)/len(mp)

# Takes list of molecules and outputs the smiles and corresponding melting point from melting_point function
collect('data/melting_point/list.txt', 'mp.txt', melting_point)
from pubchem import experimental_properties
from collect import collect
import re


def partition_coefficient(data):
    part_data = experimental_properties(data, 'Octanol/Water Partition Coefficient')

    # Octanol Water Partition Coefficient
    owpc = []

    for point in part_data['Information']:
        if 'est' not in str(point['Value']):
            owpc.append(float(re.findall(r'[-+]?[0-9]*\.?[0-9]+', str(point['Value']))[0]))

    return round(sum(owpc)/len(owpc), 2)

# Takes list of molecules and outputs the smiles and corresponding octanol/water partition coefficient from
# the partition_coefficient function
collect('list.txt', 'owpc.txt', partition_coefficient)
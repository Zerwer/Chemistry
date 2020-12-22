from pubchem import experimental_properties
from collect import collect, APIFunction
import re


class PartitionCoefficient(APIFunction):
    # APIFunction for collection LogP data
    def run(self, data):
        part_data = experimental_properties(data, 'Octanol/Water Partition Coefficient')

        # Octanol Water Partition Coefficient
        owpc = []

        for point in part_data['Information']:
            if 'est' not in str(point['Value'] and
                                'ph' not in str(point['Value'])):
                logP = float(re.findall(r'[-+]?[0-9]*\.?[0-9]+',
                                        str(point['Value']))[0])
                if logP < 500:
                    owpc.append(logP)

        return round(sum(owpc)/len(owpc), 2)


# Takes list of molecules and outputs the smiles and corresponding
#   octanol/water partition coefficient from the partition_coefficient function
collect('data/OctanolWaterPartitionCoefficient/list.txt',
        'owpc.txt', PartitionCoefficient)

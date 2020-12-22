# Modular system for scraping from PubChem
# API functions are functions that take PubChem Molecular JSON data and scrape
#   a property from it.
import datetime
import time
import sys
import json
import urllib
import urllib.error
from rdkit import RDLogger
from pubchem import process_request

# Stop garbage messages
RDLogger.DisableLog('rdApp.*')


# API functions must inherent this class
class APIFunction:
    def __init__(self):
        pass


# Consumes a PubChem mol CID and api_function and produces the output of api_
#   function when given the PubChem JSON data for the CID molecule
# api_function consumes the PubChem JSON data for a molecule
def api_get_data(cid, api_function):
    # Ensure api_function is a proper APIFunction
    assert issubclass(api_function, APIFunction)

    # Attempt to make the API request otherwise handle error
    try:
        # Download the molecule data and pass it to a function that return
        #   desired information
        link = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/'
        req = process_request(link + str(cid) + '/JSON')
        json_data = json.loads(req.decode())
        func = api_function()
        data = func.run(json_data)

        return data
    except KeyboardInterrupt:
        exit()
    except urllib.error.URLError:
        sys.stdout.write('\r')
        sys.stdout.write('-*- Connection timed out')
        sys.stdout.flush()
        time.sleep(5)
        api_get_data(cid, api_function)


# Consumes a list of molecules, an api_function and produces a new text file
#   with the name save where on each line there is a SMILES string and its
#   corresponding output from api_function
# api_function consumes the PubChem JSON data for a molecule
def collect(molecules, save, api_function):
    # Ensure api_function is a proper APIFunction
    assert issubclass(api_function, APIFunction)

    file = open(molecules, 'r')

    save = open(save, 'w')

    # Get time and print
    start_time = time.time()
    current_date = datetime.datetime.now()
    print('-*- Start Time: ' + str(current_date))

    # Size is relevant for loading bar
    lines = file.readlines()
    size = len(lines)

    for i, line in enumerate(lines):
        # Split the line, cid[0] and smiles[1]
        split_line = line.split('\t')

        # Takes the CID of the molecule and passes it through a function to
        #   get relevant data
        api_data = api_get_data(split_line[0], api_function)

        # Save data in format: MOL1 DATA1
        #                      MOL2 DATA2
        #                      .... .....
        # Smiles has \r that must be removed
        if api_data is not None:
            save.write(split_line[1][:-1] + ' ' + str(api_data) + '\n')

        # Simple loading bar code
        sys.stdout.write('\r')
        sys.stdout.write("[%-20s] %d%%" %
                         ('=' * int(20 * (i + 1) / size),
                          100 * ((i + 1) / size)))
        sys.stdout.write(' ' + str(i + 1) + '/' + str(size))
        sys.stdout.flush()

    # End of loading bar
    time_taken = round(time.time() - start_time)
    minutes = str(round((time_taken - (time_taken % 60)) / 60))
    seconds = str(time_taken % 60)
    print('\n-*- Took ' + minutes + ' minutes ' + seconds + ' seconds')

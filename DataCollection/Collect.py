# Takes list of molecules from PubChem Classification browser
# Format is:
# CID (\t) SMILES
import datetime, time, sys, json, urllib, urllib.error
from rdkit import RDLogger
from PubChem import process_request

# Stop garbage messages
RDLogger.DisableLog('rdApp.*')


# Takes cid and the function to sort through collected json data
def api_get_data(cid, api_function):
    # Attempt to make the API request otherwise handle error
    try:
        # Download the molecule data and pass it to a function that return desired information
        req = process_request('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/' + str(cid) + '/JSON')
        json_data = json.loads(req.decode())
        data = api_function(json_data)

        # This must be true unless the wrong data is used and there is no associated value
        if data:
            return data
        else:
            print('Not valid molecule: ' + str(cid))
            exit()

    except Exception as e:
        if type(e) == KeyboardInterrupt:
            exit()
        # Retry after five seconds if connection is interrupted
        elif type(e) == urllib.error.URLError:
            sys.stdout.write('\r')
            sys.stdout.write('-*- Connection timed out')
            sys.stdout.flush()
            time.sleep(5)
            api_get_data(cid, api_function)


# Function is called by function for specific data, removes redundant code
def collect(molecules, save, api_function):
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

        # Takes the CID of the molecule and passes it through a function to get relevant data
        api_data = api_get_data(split_line[0], api_function)

        # Save data in format: MOL1 DATA1
        #                      MOL2 DATA2
        # Smiles has \r that must be removed
        if api_data:
            save.write(split_line[1][:-1] + ' ' + str(api_data) + '\n')

        # Simple loading bar code
        sys.stdout.write('\r')
        sys.stdout.write("[%-20s] %d%%" % ('=' * int(20 * (i + 1) / size), 100 * ((i + 1) / size)))
        sys.stdout.write(' ' + str(i + 1) + '/' + str(size))
        sys.stdout.flush()

    # End of loading bar
    time_taken = round(time.time() - start_time)
    minutes = str(round((time_taken - (time_taken % 60)) / 60))
    seconds = str(time_taken % 60)
    print('\n-*- Took ' + minutes + ' minutes ' + seconds + ' seconds')

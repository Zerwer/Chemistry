from PubChem import *
import datetime
import urllib
import sys
from rdkit.Chem.MolStandardize.rdMolStandardize import *
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

file = open('data/boiling_point/ia_data.txt', 'r')

save = open('better.txt', 'w')

start_time = time.time()
current_date = datetime.datetime.now()
print('-*- Start Time: ' + str(current_date))

def attempt_boiling(molecule):
    try:
        pug_json = pug_smiles_json(molecule)
        cid = get_cid(pug_json)
        req = process_request(website_short + '/pug_view/data/compound/' + str(cid) + '/JSON')
        json_data = json.loads(req.decode())
        bp = boiling_point(json_data)
        if bp:
            save.write(molecule + ' ' + str(round(bp, 1)) + '\n')
    except Exception as e:
        if type(e) == KeyboardInterrupt:
            exit()
        elif type(e) == urllib.error.URLError:
            sys.stdout.write('\r')
            sys.stdout.write('-*- Connection timed out')
            sys.stdout.flush()
            time.sleep(5)
            attempt_boiling(molecule)

lines = file.readlines()

size = len(lines)

for i, line in enumerate(lines):
    attempt_boiling(line.split(' ')[0])

    sys.stdout.write('\r')
    sys.stdout.write("[%-20s] %d%%" % ('=' * int(20*((i+1))/size), 100*((i+1)/size)))
    sys.stdout.write(' ' + str(i+1) + '/' + str(size))
    sys.stdout.flush()

time_taken = round(time.time() - start_time)
minutes = str(round((time_taken - (time_taken % 60))/60))
seconds = str(time_taken % 60)

print('\n-*- Took ' + minutes + ' minutes ' + seconds + ' seconds')
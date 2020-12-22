# General functions for interacting with the PubChem API
import time
import urllib.request
import json

max_requests = 5  # Maximum requests per second
website_short = 'https://pubchem.ncbi.nlm.nih.gov/rest/'

last_request = int(round(time.time() * 1000))


# Function prevents overload from server base on max_requests
def process_request(url):
    global last_request

    current_time = int(round(time.time() * 1000))

    if current_time - last_request < 1000 / max_requests:
        time.sleep(((1000 / max_requests) -
                    (last_request - current_time)) / 1000)

    last_request = int(round(time.time() * 1000))

    return urllib.request.urlopen(urllib.request.Request(url)).read()


def pug_rest_json(identifier):
    if isinstance(identifier, int):
        req = process_request(website_short + 'pug/compound/cid/' +
                              str(identifier) + '/JSON')
    else:
        req = process_request(website_short + 'pug/compound/name/' +
                              identifier + '/JSON')

    return json.loads(req.decode())


def pug_smiles_json(smiles):
    req = process_request(website_short + 'pug/compound/smiles/' +
                          smiles + '/JSON')
    return json.loads(req.decode())


# Produces the SMILE of a molecule from JSON data
def get_smiles(json_data):
    for x in json_data['PC_Compounds'][0]['props']:
        if x['urn']['label'] == 'SMILES' and x['urn']['name'] == 'Canonical':
            return x['value']['sval']

    return ''


# Returns int CID from json data retrieved from name
def get_cid(json_data):
    return int(json_data['PC_Compounds'][0]['id']['id']['cid'])


# Takes JSON data of compound and returns data for specific
#   experimental property
def experimental_properties(data, experiment_property):
    for x in data['Record']['Section']:
        if x['TOCHeading'] == 'Chemical and Physical Properties':
            for y in x['Section']:
                if y['TOCHeading'] == 'Experimental Properties':
                    for z in y['Section']:
                        if z['TOCHeading'] == experiment_property:
                            return z
    return ''


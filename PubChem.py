import time
import urllib.request
import re
import json

max_requests = 5  # Maximum requests per second
website_short = 'https://pubchem.ncbi.nlm.nih.gov/rest/'

last_request = int(round(time.time() * 1000))

# Function prevents overload from server base on max_requests
def process_request(url):
    global last_request

    current_time = int(round(time.time() * 1000))

    if current_time - last_request < 1000 / max_requests:
        time.sleep(((1000 / max_requests) - (last_request - current_time)) / 1000)

    last_request = int(round(time.time() * 1000))

    return urllib.request.urlopen(urllib.request.Request(url)).read()


def pug_rest_json(identifier):
    if isinstance(identifier, int):
        req = process_request(website_short + 'pug/compound/cid/' + str(identifier) + '/JSON')
    else:
        req = process_request(website_short + 'pug/compound/name/' + identifier + '/JSON')

    return json.loads(req.decode())


def pug_smiles_json(smiles):
    req = process_request(website_short + 'pug/compound/smiles/' + smiles + '/JSON')
    return json.loads(req.decode())


def get_smiles(json_data):
    for x in json_data['PC_Compounds'][0]['props']:
        if x['urn']['label'] == 'SMILES' and x['urn']['name'] == 'Canonical':
            return x['value']['sval']

    return ''


# Returns int CID from json data retrieved from name
def get_cid(json_data):
    return int(json_data['PC_Compounds'][0]['id']['id']['cid'])


# Takes JSON data of compound and returns data for specific experimental property
def experimental_properties(data, experiment_property):
    for x in data['Record']['Section']:
        if x['TOCHeading'] == 'Chemical and Physical Properties':
            for y in x['Section']:
                if y['TOCHeading'] == 'Experimental Properties':
                    for z in y['Section']:
                        if z['TOCHeading'] == experiment_property:
                            return str(z)
    return ''

# -?\d+(?:\.\d+)?\s*°\s*C
def boiling_point(data):
    boiling_data = experimental_properties(data, 'Boiling Point')

    bp = []

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
                bp.append(float(re.sub('[ °C]', '', re.findall(r'-?\d+(?:\.\d+)?\s*°\s*C', line)[0])))

    if len(bp) == 0:
        return None
    else:
        return sum(bp)/len(bp)

    # temp_sum = 0.0
    # count = 0.0
    # for temp in re.findall(temp_reg, boiling_data):
    #     temp_sum += float(temp.strip(' °C'))
    #     count += 1
    #
    # if count:
    #     return round(temp_sum / count, 1)
    # else:
    #     return '?'


def get_info(identifier):
    pug_data = pug_rest_json(identifier.replace(' ', '%20'))
    cid = get_cid(pug_data)
    req = process_request(website_short + '/pug_view/data/compound/' + str(cid) + '/JSON')
    json_data = json.loads(req.decode())
    bp = boiling_point(json_data)
    smile = get_smiles(pug_data)

    print(str(identifier) + ' ' + str(bp) + ' ' + smile)

#req = process_request(website_short + '/pug_view/data/compound/' + str(123382) + '/JSON')
#req = process_request(website_short + '/pug_view/data/compound/' + str(136332) + '/JSON')
#req = process_request(website_short + '/pug_view/data/compound/' + str(15797) + '/JSON')
# req = process_request(website_short + '/pug_view/data/compound/' + str(2995) + '/JSON')
# json_data = json.loads(req.decode())
# print(boiling_point(json_data))
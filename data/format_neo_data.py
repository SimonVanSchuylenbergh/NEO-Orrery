import json
from pathlib import Path
from tqdm import tqdm


neo_names = {}
for line in open(Path(__file__).parent / 'allneo.lst', 'r').readlines():
    number = line[:8].strip()
    name = line[8:].replace('\n', '').strip()
    if number == '':
        break
    else:
        neo_names[number] = name

'''print(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][:16].strip())
print(float(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][16:32]))     # epoch
print(float(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][32:57]))     # a
print(float(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][57:82]))     # e
print(float(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][82:107]))    # i
print(float(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][107:132]))   # long node
print(float(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][132:157]))   # arg peric
print(float(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][157:181]))   # mean anomaly
print(float(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][181:187]))   # absolute magnitude
print(float(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][187:193]))   # slope param
print(bool(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6][193:]))       # non-grav param

print('---')

print(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][:27].strip())  # name
print(int(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][28:34]))   # diameter
print('*' in open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][35:44]) # diameter based on magnitude
print(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][45:63].strip()) # Impact date/time
print(float(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][64:74])) # IP max
print(float(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][75:83])) # PS max
print(int(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][84:88])) # TS
print(float(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][89:99])) # vel km/s
print(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][100:111].strip()) # years
print(float(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][112:122])) # IP cum
print(float(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4][123:131])) # PS cum'''


neo_orbital_data = {}
for line in tqdm(open(Path(__file__).parent / 'neo_kc.cat', 'r').readlines()[6:]):
    neo_orbital_data[
        neo_names[line[:16].strip()] if line[:16].strip() in neo_names else line[:16].strip()
    ] = {
        'Epoch (MJD)': float(line[16:32]),
        'a': float(line[32:57]),
        'e': float(line[57:82]),
        'inc': float(line[82:107]),
        'node': float(line[107:132]),
        'peri': float(line[132:157]),
        'ma': float(line[157:181]),
        #'slope param': float(line[181:187]),
        #'non-grav param': bool(line[193:])
    }

with open(Path(__file__).parent / 'neo_data.json', 'w') as output_file:
    output_file.write('[\n')

    for line in tqdm(open(Path(__file__).parent / 'esa_risk_list_20241004_2210.txt', 'r').readlines()[4:]):
        name = line[:27].strip()
        if not name in neo_orbital_data:
            print(name + ' not in orbital data, continuing')
            continue

        output_file.write('    {\n')
        output_file.write(f'        "name": "{name}",\n')
        output_file.write(f'        "epoch": {neo_orbital_data[name]['Epoch (MJD)']},\n')
        output_file.write(f'        "a": {neo_orbital_data[name]['a']},\n')
        output_file.write(f'        "e": {neo_orbital_data[name]['e']},\n')
        output_file.write(f'        "inc": {neo_orbital_data[name]['inc']},\n')
        output_file.write(f'        "node": {neo_orbital_data[name]['node']},\n')
        output_file.write(f'        "peri": {neo_orbital_data[name]['peri']},\n')
        output_file.write(f'        "ma": {neo_orbital_data[name]['ma']},\n')
        output_file.write(f'        "diameter": {float(line[28:34])},\n')
        output_file.write(f'        "diameter_based_on_magnitude": {str('*' in line[35:44]).lower()},\n')
        output_file.write(f'        "impact": "{line[45:63].strip()}",\n')
        output_file.write(f'        "IP max": {float(line[64:74])},\n')
        output_file.write(f'        "PS max": {float(line[75:83])},\n')
        output_file.write(f'        "TS": {int(line[84:88])},\n')
        output_file.write(f'        "vel": {float(line[89:99])},\n')
        output_file.write(f'        "years": "{line[100:111].strip()}",\n')
        output_file.write(f'        "IP cum": {float(line[112:122])},\n')
        output_file.write(f'        "PS cum": {float(line[123:131])}\n')
        output_file.write('    },\n')
    output_file.write(']\n')



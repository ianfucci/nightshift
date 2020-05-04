import argparse
from collections import defaultdict

import requests
import matplotlib.pyplot as plt

# CLI flags with descriptions
parser = argparse.ArgumentParser(description='Pull chemical shifts from BMRB and plot peak positions of amides or methyls. Currently, this only gets the data from the first entity in the entry. Not sure why...')
parser.add_argument('entry', help='BMRB accession number')
parser.add_argument('--methyl', action='store_true', help='Plots methyl shifts, not using this flag plots amide.')
parser.add_argument('-r', '--residues', help='A single string (i.e. ILV) of 1-letter amino acid codes. Will only plot these residues.')

prochiral = parser.add_mutually_exclusive_group()
prochiral.add_argument('--proR', action='store_true', help='Show only proR resonances for Leu and Val.')
prochiral.add_argument('--proS', action='store_true', help='Show only proS resonances for Leu and Val.')
args = parser.parse_args()

residue_map = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', 'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE',
               'L':'LEU', 'M':'MET', 'K':'LYS', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
residues = [residue_map[r.upper()] for r in args.residues] if args.residues else residue_map.values()

# Methyl atom names for ILVMAT
methyl_atoms = {'ILE':[('CD1','HD11')], 'LEU':[('CD1','HD11'),('CD2','HD21')], 'VAL':[('CG1','HG11'),('CG2','HG21')],
                'MET':[('CE','HE1')], 'ALA':[('CB','HB1')], 'THR':[('CG2','HG21')]}
if args.proR:
    del methyl_atoms['LEU'][1]
    del methyl_atoms['VAL'][1]
if args.proS:
    del methyl_atoms['LEU'][0]
    del methyl_atoms['VAL'][0]
amide_atoms = ('N','H')

shift_url = f'http://webapi.bmrb.wisc.edu/v2/entry/{args.entry}?saveframe_category=assigned_chemical_shifts'

# Use BMRB API to get entry number's assigned chemical shifts
with requests.get(shift_url) as r:
    bmrb_dict = r.json()

# Handle mistyped entry codes
if list(bmrb_dict.keys())[0] == 'error':
    raise KeyError(bmrb_dict['error'])

# Find where assignments are stored in BMRB file
for entry in bmrb_dict[args.entry]['assigned_chemical_shifts'][0]['loops']:
    if entry.get('category') and entry['category'] == '_Atom_chem_shift':
        bmrb_assignments = entry['data']

# Only use residue number, residue name, atom name and chemical shift
trimmed_entry = [ba[5:8] + [float(ba[10])] for ba in bmrb_assignments]
if residues:
    trimmed_entry = [te for te in trimmed_entry if te[1] in residues]

# Data structure that allows grouping of residues by atom name
# i.e. for methyls: {LeuX : {CD1: 1H_shift, 13C_shift, CD2: 1H_shift, 13C_shift}}
resonances = defaultdict(lambda: defaultdict(list))
for num,res,atom,shift in trimmed_entry:
    if not args.methyl:
        if atom in amide_atoms:
            resonances[res+num]['NH'].append(shift)
    if args.methyl:
        pairs = methyl_atoms.get(res)
        if pairs is None:
            continue
        for p in pairs:
            if atom in p:
                resonances[res+num][p[0]].append(shift)

# Plot points on "simulated spectrum"
fig = plt.figure(figsize=(12,12))
plt.xlabel('1H (ppm)')
f1 = '13C (ppm)' if args.methyl else '15N (ppm)'
plt.ylabel(f1)

for res, atom in resonances.items():
    for name, point in atom.items():
        try:
            # Make heteroatom the y axis
            if point[0] > point[1]:
                point = point[::-1]
            plt.scatter(*point)
            # Label each point on the spectrum
            plt.annotate(res+name, point)
        except (IndexError,TypeError):
            # Skip if only one shift is given
            continue

# Invert axes
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

plt.show()
plt.close()

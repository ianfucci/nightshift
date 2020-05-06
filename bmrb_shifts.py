import argparse
import warnings
from collections import defaultdict

import requests
import matplotlib.pyplot as plt

# CLI flags with descriptions
parser = argparse.ArgumentParser(description='Pull assigned chemical shifts from BMRB entry and plot peak positions on a simulated 2D "spectrum".')
parser.add_argument('entry', help='BMRB accession number')
parser.add_argument('-r', '--residues', help='Residues filter: a single string (i.e. ILV) of 1-letter amino acid codes. Will only plot these residues.')

# Flags to set which correlation will be depicted
correlation = parser.add_mutually_exclusive_group()
correlation.add_argument('--amide', action='store_true', help='Plot amide shifts.')
correlation.add_argument('--methyl', nargs='*', choices=['proR', 'proS'], help='Plot methyl shifts, allows use of proR and proS flags for specific labeling schemes. By default plots both.')
correlation.add_argument('--custom', nargs=2, help='Plot correlation between any two atoms (i.e. N CO).')
args = parser.parse_args()
args.custom = [a.upper() for a in args.custom]

residue_map = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', 'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE',
               'L':'LEU', 'M':'MET', 'K':'LYS', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
# Methyl atom names for MILVAT
methyl_atoms = {'ILE':[('CD1','HD11')], 'LEU':[('CD1','HD11'),('CD2','HD21')], 'VAL':[('CG1','HG11'),('CG2','HG21')],
                'MET':[('CE','HE1')], 'ALA':[('CB','HB1')], 'THR':[('CG2','HG21')]}

# Check residue filter
if args.residues:
    # Warn for incorrect 1-letter codes
    bad_codes = set(args.residues.upper()) - residue_map.keys()
    if bad_codes:
        bad_code_string = ','.join(bad_codes)
        warnings.warn(f'{bad_code_string} not valid 1-letter code(s)')
    # Remove bad codes and ignore
    residues = [residue_map.get(r.upper()) for r in args.residues if r not in bad_codes]
    # Warn if non MILVAT residues are used with the --methyl flag
    non_milvat = set(residues) - methyl_atoms.keys()
    if args.methyl is not None and non_milvat:
        res_string = ','.join(non_milvat)
        warnings.warn(f'residues other than MILVAT: ({res_string}) were ignored when plotting methyl spectra')
else:
    # If no residue filter is specified use all residues
    residues = residue_map.values()

# This is the same as the methyl argument and this makes it easier to handle
# Cmethyl and Hmethyl functionality not yet implemented
if 'CMETHYL' in args.custom and 'HMETHYL' in args.custom:
    args.methyl = []
    args.custom = None
# Make selector for Cmethyl and Hmethyl special arguments for custom
elif 'CMETHYL' in args.custom:
    c_selector = {}
    for res, pairs in methyl_atoms.items():
        c_selector[res] = [p[0] for p in pairs]
elif 'HMETHYL' in args.custom:
    h_selector = {}
    for res, pairs in methyl_atoms.items():
        h_selector[res] = [p[1] for p in pairs]

if args.amide:
    # Make a dictionary similar to methyl_atoms
    amide_atoms = dict(zip(residues, [[('N','H')]]*len(residues)))
    # Selector is a function that returns pairs of atoms given a residue
    selector = amide_atoms
# args.methyl will be an empty list for default or a list ['proR'] or ['proS']
# will be None for --amide or --custom
elif args.methyl is not None:
    for res, pairs in methyl_atoms.items():
        if res in {'LEU','VAL'} and args.methyl:
            # Filter out prochiral atoms when given proR/proS flags
            if args.methyl[0] == 'proR':
                pairs = [p for p in pairs if not p[0].endswith('2')]
                methyl_atoms[res] = pairs
            elif args.methyl[0] == 'proS':
                pairs = [p for p in pairs if not p[0].endswith('1')]
                methyl_atoms[res] = pairs
    # Selector is a function that returns pairs of atoms given a residue
    selector = methyl_atoms

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

# Only use residue number, residue name, atom name and chemical shift and filter residues
assignments = defaultdict(dict)
for ba in bmrb_assignments:
    if ba[6] in residues:
        # format {RES<num>: {atom1: shift, atom2: shift}}
        assignments[ba[6]+ba[5]][ba[7]] = float(ba[10])
        print(ba[6],ba[7])

# Data structure containing selected atoms' peaks
resonances = defaultdict(dict)
for name, atom_dict in assignments.items():
    res = name[:3]
    # Custom atom names run into a lot of different cases
    if args.custom is not None:
        simple = {'H','N','C','CA','CB'}
        # Simple case, just index atom dict based on atom name
        if args.custom[0] in simple and args.custom[1] in simple:
            resonances[name]['_'.join(args.custom)] = (atom_dict.get(args.custom[1]), atom_dict.get(args.custom[0]))
        # Only the first argument is the simple case
        elif args.custom[0] in simple and args.custom[1] not in simple:
            second = atom_dict.get(args.custom[0])
            for atom, shift in atom_dict.items():
                if atom.startswith(args.custom[1]):
                    resonances[name][args.custom[0]+'_'+atom] = (shift, second)
        # Only the second argument is the simple case
        elif args.custom[1] in simple and args.custom[0] not in simple:
            first = atom_dict.get(args.custom[1])
            for atom, shift in atom_dict.items():
                if atom.startswith(args.custom[0]):
                    resonances[name][atom+'_'+args.custom[1]] = (first, shift)
        # Both need to match the first letters of the atom ID
        else:
            # Iterate through twice to capture all matches
            for atom, shift in atom_dict.items():
                if atom.startswith(args.custom[0]):
                    for atom2, shift2 in atom_dict.items():
                        if atom2.startswith(args.custom[1]):
                            resonances[name][atom+'_'+atom2] = (shift2, shift)
    else:
        # amide and methyl are special cases, used predefined selectors
        pairs = selector.get(res)
        if pairs is not None:
            for p in pairs:
                # indices are swapped for plotting purposes, prevents doing it later
                resonances[name]['_'.join(p)] = (atom_dict.get(p[1]), atom_dict.get(p[0]))

# Plot points on "simulated spectrum"
fig = plt.figure()

# Axes labels hardcoded for amide and methyl
if args.methyl is not None:
    plt.ylabel('13C (ppm)')
    plt.xlabel('1H (ppm)')
elif args.amide: 
    plt.ylabel('15N (ppm)')
    plt.xlabel('1H (ppm)')
# Axes labels for custom are the atom IDs
else:
    if args.custom[0].startswith('H'):
        nucleus1 = '1H'
    elif args.custom[0].startswith('C'):
        nucleus1 = '13C'
    elif args.custom[0].startswith('N'):
        nucleus1 = '15N'
    if args.custom[1].startswith('H'):
        nucleus2 = '1H'
    elif args.custom[1].startswith('C'):
        nucleus2 = '13C'
    elif args.custom[1].startswith('N'):
        nucleus2 = '15N'
    plt.ylabel(f'{args.custom[0]} - {nucleus1} (ppm)')
    plt.xlabel(f'{args.custom[1]} - {nucleus2} (ppm)')

# Palette of distinct colors for each residue
residue_colors = {'ALA':'maroon', 'ARG':'red', 'ASN':'pink', 'ASP':'brown', 'CYS':'orange', 'GLU':'coral', 'GLN':'olive', 'GLY':'magenta', 'HIS':'khaki', 'ILE':'purple',
                  'LEU':'green', 'MET':'navy', 'LYS':'blue', 'PHE':'lime', 'PRO':'lightgreen', 'SER':'aquamarine', 'THR':'cyan', 'TRP':'black', 'TYR':'grey', 'VAL':'yellow'}

for name, pairs in resonances.items():
    res = name[:3]
    for atoms, peak in pairs.items():
        if None not in peak:
            # Plot peak centers colored by residue
            plt.scatter(*peak, c=residue_colors[res])
            # Label each peak, add atom ID to prochiral atoms
            if (args.methyl is not None and res in {'LEU', 'VAL'}) or args.custom is not None:
                peak_label = name + atoms.split('_')[0]
                plt.annotate(peak_label, peak)
            else:
                plt.annotate(name, peak)

# Invert axes
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

plt.show()
plt.close()

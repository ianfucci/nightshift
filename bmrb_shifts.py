#!/usr/bin/env python

'''This is the main script to produce simulated 2D spectra for protein entries from the BMRB
'''

import argparse
import itertools
import warnings
from collections import defaultdict

import requests
import matplotlib.pyplot as plt

# CLI flags with descriptions
parser = argparse.ArgumentParser(description='Pull assigned chemical shifts from BMRB entry and plot peak positions on a simulated 2D "spectrum".')
parser.add_argument('entry', help='BMRB accession number')
parser.add_argument('-r', '--residues', help='Residues filter: a single string (i.e. ILV) of 1-letter amino acid codes. Will only plot these residues.')
parser.add_argument('--offset', default='0', help='Index offset. BMRB indices starts at 1 which may not reflect the true residue number')
parser.add_argument('--csv', help='Name for output file of shifts in CSV format for later plotting.')
parser.add_argument('--nolabels', action='store_true', help='Plot spectrum without peak labels.')
parser.add_argument('--showlegend', action='store_true', help='Show a legend on the spectrum with residue colors.')
parser.add_argument('-o', '--output', help='Save spectrum to file and do not preview plot with matplotlib.')

# Flags to set which correlation will be depicted
correlation = parser.add_mutually_exclusive_group()
correlation.add_argument('--amide', action='store_true', help='Plot amide shifts.')
correlation.add_argument('--methyl', nargs='*', choices=['proR', 'proS'], help='Plot methyl shifts, allows use of proR and proS flags for specific labeling schemes. By default plots both.')
correlation.add_argument('--custom', nargs=2, help='Plot correlation between any two atoms (i.e. N CO).')
args = parser.parse_args()

residue_map = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU',
               'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'M':'MET',
               'K':'LYS', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP',
               'Y':'TYR', 'V':'VAL'}
# Methyl atom names for MILVAT
methyl_atoms = {'ILE':[('CD1', 'HD11')], 'LEU':[('CD1', 'HD11'), ('CD2', 'HD21')],
                'VAL':[('CG1', 'HG11'), ('CG2', 'HG21')], 'MET':[('CE', 'HE1')],
                'ALA':[('CB', 'HB1')], 'THR':[('CG2', 'HG21')]}
simple_atoms = {'H', 'HA', 'N', 'C', 'CA', 'CB'}

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
        warnings.warn(f'residues other than MILVAT: ({res_string}) are ignored when plotting')
else:
    # If no residue filter is specified use all residues
    residues = residue_map.values()

# Perform correlation specific setup
minus_one = False
if args.custom is not None:
    args.custom = [a.upper() for a in args.custom]
    # This is the same as the methyl argument and this makes it easier to handle
    if 'CMETHYL' in args.custom and 'HMETHYL' in args.custom:
        args.methyl = []
        args.custom = None
    for i, custom_atom in enumerate(args.custom):
        if custom_atom.endswith('-1'):
            minus_one = True
            # Save which atom choice has the (i-1)
            args.custom[i] = custom_atom.split('-')[0]
            prev_index = i
if args.amide:
    # Make a dictionary similar to methyl_atoms
    selector = dict(zip(residues, [[('N', 'H')]]*len(residues)))
# args.methyl will be an empty list for default or a list ['proR'] or ['proS']
# will be None for --amide or --custom
if args.methyl is not None:
    for res, pairs in methyl_atoms.items():
        if res in {'LEU', 'VAL'} and args.methyl:
            # Filter out prochiral atoms when given proR/proS flags
            if args.methyl[0] == 'proR':
                pairs = [p for p in pairs if not p[0].endswith('2')]
                methyl_atoms[res] = pairs
            elif args.methyl[0] == 'proS':
                pairs = [p for p in pairs if not p[0].endswith('1')]
                methyl_atoms[res] = pairs
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

# Parse BMRB assignments file storing relevant fields
# Only use residue number, residue name, atom name and chemical shift
assignments = defaultdict(dict)
encountered_atoms = defaultdict(list)
for ba in bmrb_assignments:
    # format for assignments dict {RES<num>: {atom1: shift, atom2: shift...}}
    assignments[ba[6]+ba[5]][ba[7]] = float(ba[10])
    if encountered_atoms.get(ba[6]) is None or ba[7] not in encountered_atoms[ba[6]]:
        encountered_atoms[ba[6]].append(ba[7])

# Generate selectors for custom correlations from encountered atom names
if args.custom is not None:
    selector_dicts = []
    for custom_atom in args.custom:
        # Cmethyl and Hmethyl already defined by modifying methyl_atoms
        if custom_atom == 'CMETHYL':
            c_selector = {}
            for res, pairs in methyl_atoms.items():
                c_selector[res] = [p[0] for p in pairs]
            selector_dicts.append(c_selector)
        elif custom_atom == 'HMETHYL':
            h_selector = {}
            for res, pairs in methyl_atoms.items():
                h_selector[res] = [p[1] for p in pairs]
            selector_dicts.append(h_selector)
        # If the atom is not a simple case match against the first letters of the atom IDs
        elif custom_atom not in simple_atoms:
            custom_selector = {}
            for res, atoms in encountered_atoms.items():
                custom_selector[res] = [a for a in atoms if a.startswith(custom_atom)]
            selector_dicts.append(custom_selector)
        # custom_atom is the simple case
        else:
            simple_selector = dict(zip(residues, [[(custom_atom)]]*len(residues)))
            selector_dicts.append(simple_selector)
    selector = {}
    # Does not matter which is longer, will either be filtered here or when plotting
    for res in selector_dicts[0]:
        try:
            selector[res] = list(itertools.product(selector_dicts[0][res], selector_dicts[1][res]))
        except KeyError:
            # Key was not present in both dicts
            continue

# Data structure containing selected atoms' peaks
resonances = defaultdict(dict)
# For an i-1 correlation
if minus_one:
    # Split into two iterators prev is residue (i-1) curr is residue (i)
    prev, curr = itertools.tee(assignments.items(), 2)
    curr = itertools.islice(curr, 1, None)
    lookbehind = itertools.zip_longest(prev, curr, fillvalue=None)

    curr_index = 0 if prev_index == 1 else 1
    for prev_assign, curr_assign in lookbehind:
        try:
            prev_atom_dict = prev_assign[1]
            curr_atom_dict = curr_assign[1]
            curr_res = curr_assign[0][:3]
            if int(curr_assign[0][3:]) - 1 == int(prev_assign[0][3:]):
                pairs = selector.get(curr_res)
                if curr_res in residues and pairs is not None:
                    for p in pairs:
                        if curr_index < prev_index:
                            resonances[curr_assign[0]]['_'.join(p)] = (prev_atom_dict.get(p[prev_index]), curr_atom_dict.get(p[curr_index]))
                        else:
                            resonances[curr_assign[0]]['_'.join(p)] = (curr_atom_dict.get(p[curr_index]), prev_atom_dict.get(p[prev_index]))
        except TypeError:
            # Reached None the end of the curr iterator
            continue
# For a typical spectrum
else:
    for name, atom_dict in assignments.items():
        res = name[:3]
        pairs = selector.get(res)
        # Apply residue filter
        if res in residues and pairs is not None:
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
residue_colors = {'ALA':'maroon', 'ARG':'red', 'ASN':'pink', 'ASP':'brown',
                  'CYS':'orange', 'GLU':'coral', 'GLN':'olive', 'GLY':'magenta',
                  'HIS':'khaki', 'ILE':'purple', 'LEU':'green', 'MET':'navy',
                  'LYS':'blue', 'PHE':'lime', 'PRO':'lightgreen', 'SER':'aquamarine',
                  'THR':'cyan', 'TRP':'black', 'TYR':'grey', 'VAL':'yellow'}

# Create output file
if args.csv:
    output = open(args.csv, 'w')

legend = {}
for name, pairs in resonances.items():
    res = name[:3]
    idx = int(name[3:]) + int(args.offset)
    for atoms, peak in pairs.items():
        if None not in peak:
            # Plot peak centers colored by residue
            handle = plt.scatter(*peak, c=residue_colors[res])
            legend[res] = handle
            # Label each peak, add atom ID to prochiral atoms
            if (args.methyl is not None and res in {'LEU', 'VAL'}) or args.custom is not None:
                peak_label = res + str(idx) + atoms.split('_')[0]
                if not args.nolabels:
                    plt.annotate(peak_label, peak)
                if args.csv:
                    output.write(f'{peak_label},' + ','.join(str(p) for p in peak) + '\n')
            else:
                if not args.nolabels:
                    plt.annotate(f'{res}{idx}', peak)
                if args.csv:
                    output.write(f'{res}{idx},' + ','.join(str(p) for p in peak) + '\n')

if args.nolabels or args.showlegend:
    sort_legend = {k: v for k, v in sorted(legend.items())}
    plt.legend(sort_legend.values(), sort_legend.keys())
    
# Invert axes
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

if args.csv:
    output.close()

if not args.output:
    plt.show()
else:
    plt.savefig(args.output, dpi=300)
plt.close()

import argparse
from collections import namedtuple
import csv
import logging
import re
import sys

import matplotlib.pyplot as plt

from nightshift import bmrb
from nightshift import constants
from nightshift import plotting
from nightshift.selector import Selector, AmideSelector, MethylSelector, AdvancedSelector

def run_cli() -> None:
    parser = argparse.ArgumentParser(description='Tools for plotting simulated NMR spectra from assigned chemical shifts in the BMRB.')
    subparsers = parser.add_subparsers()

    # Subcommand `get` is the primary function of the script
    get_parser = subparsers.add_parser('get', description='Pull assigned chemical shifts from BMRB entry and plot peak positions on a simulated spectrum.')
    get_parser.add_argument('entry', help='BMRB entry number. If not knonw, use `nightshift search`')
    get_parser.add_argument('-r', '--residues', help='Residues filter: a single string (i.e. ILV) of 1-letter amino acid codes. Will only plot these residues.')
    get_parser.add_argument('--offset', default=0, type=int, help='Index offset. BMRB indices starts at 1 which may not reflect the true residue number.')
    get_parser.add_argument('--segment', nargs=2, type=int, help='Indices representing the start and stop of the sequence you want to show. Does not account for offset Format: start stop')
    get_parser.add_argument('--csv', help='Name for output file of shifts in CSV format for later plotting.')
    get_parser.add_argument('--label', type=int, help='Choose which residue is used for the label in custom correlations: 1, 2 or 3. Default is 1.')
    get_parser.add_argument('--project', default=3, type=int, help='Choose which dimension to project on for 3D spectra, 1, 2 or 3.')
    get_parser.add_argument('--slices', default=16, type=int, help='Number of slices/points in 3rd dimension in 3D spectra. Default is 16')
    get_parser.add_argument('--nolabels', action='store_true', help='Plot spectrum without peak labels.')
    get_parser.add_argument('--showlegend', action='store_true', help='Show a legend on the spectrum with residue colors.')
    get_parser.add_argument('-o', '--output', help='Save spectrum to file and do not plot in interactive mode with matplotlib.')
    
    # Flags to set which correlation will be depicted
    correlation = get_parser.add_mutually_exclusive_group(required=True)
    correlation.add_argument('--amide', action='store_true', help='Plot amide shifts.')
    correlation.add_argument('--methyl', action='store_true', help='Plot methyl shifts, allows use of --proR and --proS flags for specific labeling schemes. By default plots both.')
    correlation.add_argument('--custom', nargs='*', help='Plot correlation between any two or three atoms (i.e. N C). Including +<number> or -<number> after an atom name allows for inter-residue correlations (i.e. H N C-1) Also supports HMETHYL and CMETHYL and proR/proS.')

    atom_filters = get_parser.add_mutually_exclusive_group()
    atom_filters.add_argument('--sidechains', action='store_true', help='Choose to plot side chains of Asn and Gln. Works with --amide')
    atom_filters.add_argument('--proR', action='store_true', help='Choose to plot only proR peaks for Leu and Val. Works with --methyl and --custom')
    atom_filters.add_argument('--proS', action='store_true', help='Choose to plot only proS peaks for Leu and Val. Works with --methyl and --custom')
    get_parser.set_defaults(func=get)
    
    # Allows spectra saved to CSV files to be replotted or superimposed
    open_parser = subparsers.add_parser('open', description='Plot one or more spectra from CSV output files saved by `get`. Allows spectra saved to CSV files to be replotted or superimposed.')
    open_parser.add_argument('input', nargs='*', help='Name of input CSV files to overlay.')
    open_parser.add_argument('--csv', help='Name for output file of shifts in CSV format for later plotting.')
    open_parser.add_argument('--nolabels', action='store_true', help='Plot spectrum without peak labels.')
    open_parser.add_argument('--showlegend', action='store_true', help='Show a legend on the spectrum with residue colors.')
    open_parser.add_argument('-o', '--output', help='Save spectrum to file and do not preview plot with matplotlib.')
    open_parser.set_defaults(func=from_file)

    # Search BMRB for an entry ID using text or a PDB ID
    search_parser = subparsers.add_parser('search', description='Search the BMRB for an entry number for use with `get`. Can use simple search terms including a PDB code.')
    search_parser.add_argument('terms', nargs='*', help='Search terms for BMRB database.')
    search_parser.set_defaults(func=search)

    # Open BMRB webpage for an entry number
    website_parser = subparsers.add_parser('website', description='Open BMRB webpage for an entry number')
    website_parser.add_argument('entry', help='BMRB entry number. If not knonw, use `nightshift search`')
    website_parser.set_defaults(func=website)

    # Parse arguments and run assigned function
    args = parser.parse_args()
    args.func(args)

def get(args):
    # Check residue filter
    if args.residues:
        # Warn for incorrect 1-letter codes
        bad_codes = set(args.residues.upper()) - constants.ONE_LETTER_TO_THREE_LETTER.keys()
        if bad_codes:
            bad_code_string = ','.join(bad_codes)
            logging.warn(f'{bad_code_string} not valid 1-letter code(s)')
        # Remove bad codes and ignore
        residues = [constants.ONE_LETTER_TO_THREE_LETTER.get(r.upper()) for r in args.residues if r not in bad_codes]
        # Warn if non MILVAT residues are used with the --methyl flag
        non_milvat = set(residues) - constants.METHYL_ATOMS.keys()
        if args.methyl and non_milvat:
            res_string = ','.join(non_milvat)
            logging.warn(f'residues other than MILVAT: ({res_string}) are ignored when plotting.')
    else:
        # If no residue filter is specified use all residues
        residues = list(constants.ONE_LETTER_TO_THREE_LETTER.values())

    # Get data from BMRB
    entry_data = bmrb.get_bmrb_shifts(args.entry)
    if entry_data is None:
        sys.exit(1)

    # Generate selector objects
    if args.amide:
        selector = AmideSelector(residues, args.segment, sidechains=args.sidechains)
        correlations = selector.get_correlations(entry_data)
    elif args.methyl:
        selector = MethylSelector(residues, args.segment, proR=args.proR, proS=args.proS)
        correlations = selector.get_correlations(entry_data)
    elif all(atom in constants.SIMPLE_ATOMS for atom in args.custom):
        selector = Selector(residues, args.custom, args.segment)
        correlations = selector.get_correlations(entry_data)
    else:
        # Advanced correlation setup
        args.custom = [a.upper() for a in args.custom]
        plus_minus = [0]*len(args.custom)
        atoms = []
        for i, spin in enumerate(args.custom):
            try:
                # Split at plus or minus sign
                atom, sign, index = re.split(r'(\+|\-)', spin)
                atoms.append(atom)
                plus_minus[i] = int(sign+index)
            except ValueError:
                # No plus or minus for atom
                atoms.append(spin)
                continue
        # Ensure there is an i residue, not all have +/- indices
        if 0 not in plus_minus:
            abs_min = min([abs(pm) for pm in plus_minus])
            add_min = [pm + abs_min for pm in plus_minus]
            if 0 in add_min:
                plus_minus = add_min
            else:
                plus_minus = [pm - abs_min for pm in plus_minus]
            logging.warn(f"No 'i-residue' found. Adjusted by indicies by {abs_min}.")
        selector = AdvancedSelector(residues, tuple(atoms), args.segment, plus_minus=plus_minus, proR=args.proR, proS=args.proS)
    
        if args.label is not None:
            args.label -= 1
        correlations = selector.get_correlations(entry_data, label=args.label)
    
    if args.csv is not None:
        with open(args.csv, 'w', newline='') as csvfile:
            fieldnames = ['label'] + list(selector.atoms)
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, quotechar="'", quoting=csv.QUOTE_NONNUMERIC)
            writer.writeheader()
            for sequence_number, residue_type, chemical_shifts in correlations:
                row = {atom : chemical_shifts[i] for atom in enumerate(selector.atoms)}
                row['label'] = f'{residue_type}{sequence_number + args.offset}'
                writer.writerow(row)

    fig, ax = plt.subplots()
    # Plot correlations
    if len(selector.atoms) == 2:
        plotting.plot2D(ax, correlations, selector.atoms, nolabels=args.nolabels, showlegend=args.showlegend, offset=args.offset)

    elif len(selector.atoms) == 3:
        # Have to return tracker, otherwise gets GC'd and scrolling doesn't work.
        tracker = plotting.plot3D(ax, correlations, selector.atoms, nolabels=args.nolabels, showlegend=args.showlegend, offset=args.offset,
        project=args.project-1, slices=args.slices)

    # Interactive matplotlib window opened if not saving
    if not args.output:
        plt.show()
    else:
        plt.savefig(args.output, dpi=300)
    plt.close()


def from_file(args):
    records = []
    for in_file in args.input:
        with open(in_file, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile, quotechar="'", quoting=csv.QUOTE_NONNUMERIC)
            NightshiftRecord = namedtuple('NightshiftRecord', reader.fieldnames)
            for row in reader:
                records.append(NightshiftRecord(*row.values()))
    if all(len(record) == 3 for record in records):
        # 2 + 1 for label
        plotting.plot2D(records, nolabels=args.nolabels, showlegend=args.showlegend, segment=(310,319), offset=0)
    elif all(len(record) == 4 for record in records):
        # 3 + 1 for label
        pass
    else:
        pass

    # Interactive matplotlib window opened if not saving
    if not args.output:
        plt.show()
    else:
        plt.savefig(args.output, dpi=300)
    plt.close()

def search(args: argparse.Namespace) -> None:
    if args.terms:
        terms = ' '.join(args.terms)
    else:
        terms = ''
    bmrb.bmrb_instant_search(terms)

def website(args: argparse.Namespace) -> None:
    bmrb.open_bmrb_page(args.entry)
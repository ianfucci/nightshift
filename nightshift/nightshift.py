import argparse
from collections import namedtuple
import csv
from typing import List, Tuple
import sys

import matplotlib.pyplot as plt

from nightshift import bmrb
from nightshift import constants
from nightshift import plotting
from nightshift.selector import Selector, AmideSelector, MethylSelector, AdvancedSelector

Correlations = List[Tuple[int,str,Tuple[float]]]

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
    # Get data from BMRB
    entry_data = bmrb.get_bmrb_shifts(args.entry)
    if entry_data is None:
        sys.exit(1)

    # Generate selector objects
    if args.amide:
        selector = AmideSelector(args.residues, args.segment, sidechains=args.sidechains)
        correlations = selector.get_correlations(entry_data)
    elif args.methyl:
        selector = MethylSelector(args.residues, args.segment, proR=args.proR, proS=args.proS)
        correlations = selector.get_correlations(entry_data)
    elif all(atom in constants.SIMPLE_ATOMS for atom in args.custom):
        selector = Selector(args.residues, args.custom, args.segment)
        correlations = selector.get_correlations(entry_data)
    else:
        selector = AdvancedSelector(args.residues, args.custom, args.segment, proR=args.proR, proS=args.proS)
    
        if args.label is not None:
            args.label -= 1
        correlations = selector.get_correlations(entry_data, label=args.label)
    
    if args.csv is not None:
        write_csv(args.csv, correlations, selector.atoms, args.offset)

    _, ax = plt.subplots()
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
    correlations = []
    for in_file in args.input:
        with open(in_file, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
            atoms = tuple(next(reader))[1:] # get headings
            for row in reader:
                # sequence_number, residue_type, (shift1, shift2,...)
                correlations.append([int(row[0][3:]), row[0][:3], tuple(row[1:])])

    _, ax = plt.subplots()
    if all(len(correlation[-1]) == 2 for correlation in correlations):
        plotting.plot2D(ax, correlations, atoms, nolabels=args.nolabels, showlegend=args.showlegend)
    elif all(len(correlation[-1]) == 3 for correlation in correlations):
        tracker = plotting.plot3D(ax, correlations, atoms, nolabels=args.nolabels, showlegend=args.showlegend)


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

def write_csv(filename: str, correlations: Correlations, atoms: Tuple[str], offset: int):
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = ['label'] + list(atoms)
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, quoting=csv.QUOTE_NONNUMERIC)
        writer.writeheader()
        for sequence_number, residue_type, chemical_shifts in correlations:
            row = {atom : chemical_shifts[i] for i, atom in enumerate(atoms)}
            row['label'] = f'{residue_type}{sequence_number + offset}'
            writer.writerow(row)
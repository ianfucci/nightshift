import argparse
import csv
import itertools
from typing import List, Tuple
import sys

import matplotlib.pyplot as plt

from nightshift import bmrb
from nightshift import cli
from nightshift import constants
from nightshift import plotting
from nightshift.selector import Selector, AmideSelector, MethylSelector, AdvancedSelector

Correlations = List[Tuple[int,str,Tuple[float]]]

def main() -> None:
    # get name of this module rather than __main__.py
    prog = __name__.rpartition('.')[-1]
    parser = cli.build_parser(prog)
    args = parser.parse_args()
    try:
        # run functions assigned to keywords
        args.func(args)
    except AttributeError:
        # no arguments passed
        parser.print_help()

def get(args: argparse.Namespace) -> None:
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
        plotting.plot2D(ax, correlations, selector.atoms, color=args.color, nolabels=args.nolabels, showlegend=args.showlegend, offset=args.offset)

    elif len(selector.atoms) == 3:
        # Have to return Slices3D object, otherwise gets GC'd and scrolling doesn't work.
        _ = plotting.plot3D(ax, correlations, selector.atoms, color=args.color, nolabels=args.nolabels, showlegend=args.showlegend, offset=args.offset,
        project=args.project-1, slices=args.slices)

    # Interactive matplotlib window opened if not saving
    if not args.output:
        plt.show()
    else:
        plt.savefig(args.output, dpi=300)
    plt.close()

def from_file(args: argparse.Namespace) -> None:
    _, ax = plt.subplots()
    color_cycle = itertools.cycle(args.colors)
    
    for in_file, color in zip(args.input, color_cycle):
        correlations = []
        with open(in_file, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
            atoms = tuple(next(reader))[1:] # get headings
            for row in reader:
                # sequence_number, residue_type, (shift1, shift2,...)
                correlations.append([int(row[0][3:]), row[0][:3], tuple(row[1:])])

        if all(len(correlation[-1]) == 2 for correlation in correlations):
            plotting.plot2D(ax, correlations, atoms, color=color, nolabels=args.nolabels, showlegend=args.showlegend)
        
        elif all(len(correlation[-1]) == 3 for correlation in correlations):
            # Have to return Slices3D object, otherwise gets GC'd and scrolling doesn't work.
            _ = plotting.plot3D(ax, correlations, atoms, color=color, nolabels=args.nolabels, showlegend=args.showlegend)

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

def write_csv(filename: str, correlations: Correlations, atoms: Tuple[str], offset: int) -> None:
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = ['label'] + list(atoms)
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, quoting=csv.QUOTE_NONNUMERIC)
        writer.writeheader()
        for sequence_number, residue_type, chemical_shifts in correlations:
            row = {atom : chemical_shifts[i] for i, atom in enumerate(atoms)}
            row['label'] = f'{residue_type}{sequence_number + offset}'
            writer.writerow(row)
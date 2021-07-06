import argparse
import csv
import itertools
import shutil
import sys
from typing import List, TextIO, Tuple

import matplotlib.pyplot as plt

from nightshift import bmrb
from nightshift import cli
from nightshift import constants
from nightshift import plotting
from nightshift.selector import Selector, AmideSelector, MethylSelector, AdvancedSelector, GroupSelector

Correlations = List[Tuple[int,str,Tuple[float]]]

def main() -> None:
    # get name of this module rather than __main__.py
    prog = __name__.rpartition('.')[-1]
    parser = cli.build_parser(prog)
    args = parser.parse_args()
    try:
        # check if arguments were passed and function was assigned
        getattr(args, 'func')
    except AttributeError:
        # no arguments passed
        parser.print_help()
    else:
        # run functions assigned to keywords
        args.func(args)

def get(args: argparse.Namespace) -> None:
    entry_data = bmrb.get_shifts(args.entry)

    if entry_data is None:
        sys.exit(1)

    entity_selection = 0
    if len(entry_data.shifts) > 1:
        entity_selection = choose_entity(entry_data.names)
    entity_shifts = entry_data.shifts[entity_selection]

    # Generate selector objects
    if args.amide:
        selector = AmideSelector(args.residues, args.segment, sidechains=args.sidechains)
        correlations = selector.get_correlations(entity_shifts)
    
    elif args.methyl:
        selector = MethylSelector(args.residues, args.segment, proR=args.proR, proS=args.proS)
        correlations = selector.get_correlations(entity_shifts)
    
    elif all(atom in constants.SIMPLE_ATOMS for atom in args.custom):
        # This fails for atom groups, which is intended
        selector = Selector(args.residues, args.custom, args.segment)
        correlations = selector.get_correlations(entity_shifts)
    
    elif all(isinstance(atom, tuple) for atom in args.custom):
        selector = GroupSelector(args.residues, args.custom, args.segment, proR=args.proR, proS=args.proS)
        if args.label is not None:
            args.label -= 1
        correlations = selector.get_correlations(entity_shifts, label=args.label)
    
    else:
        selector = AdvancedSelector(args.residues, args.custom, args.segment, proR=args.proR, proS=args.proS)
        if args.label is not None:
            args.label -= 1
        correlations = selector.get_correlations(entity_shifts, label=args.label)
    
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
        with in_file as csvfile:
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
    
    results = bmrb.instant_search(terms)
    if not results or results is None:
        print(f"No entries with assigned chemical shifts for '{terms}'")
        sys.exit(1)

    padding = 8
    terminal_width = shutil.get_terminal_size()[0]
    for entry_id, title in results:
        line = f"{entry_id}\t{title}"
        if len(line) + padding > terminal_width:
            line = line[:terminal_width - padding] + '...'
        print(line)

def website(args: argparse.Namespace) -> None:
    bmrb.open_website(args.entry)

def write_csv(outfile: TextIO, correlations: Correlations, atoms: Tuple[str], offset: int) -> None:
    with outfile as csvfile:
        fieldnames = ['label'] + list(atoms)
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, quoting=csv.QUOTE_NONNUMERIC)
        writer.writeheader()
        for sequence_number, residue_type, chemical_shifts in correlations:
            row = {atom : chemical_shifts[i] for i, atom in enumerate(atoms)}
            row['label'] = f'{residue_type}{sequence_number + offset}'
            writer.writerow(row)

def choose_entity(names: list) -> int:
    while True:
        try:
            print('Found more than one assigned chemical shift list, please select:')
            # Print list of names received from names_url request number matches Entity_ID tag
            for i, name in enumerate(names, start=1):
                print(f'({i}) {name}')
            entity_selection = int(input('> ')) - 1 # Entity_ID is 1-indexed, python is 0-indexed
            if 0 <= entity_selection < len(names):
                return entity_selection
            else:
                print(f'Invalid selection please choose a number <= {len(names)}.')
        except ValueError:
            print(f'Invalid selection, please input a number <= {len(names)}.')
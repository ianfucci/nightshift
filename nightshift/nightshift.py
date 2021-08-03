import argparse
import csv
from collections import defaultdict
import itertools
import logging
import shutil
import sys
from typing import List, TextIO, Tuple

import matplotlib.pyplot as plt
plt.switch_backend('Qt5Agg')

from nightshift import bmrb
from nightshift import cli
from nightshift import constants
from nightshift import plotting
from nightshift.selector import Selector, AmideSelector, MethylSelector, AdvancedSelector, GroupSelector

Correlations = List[Tuple[int,str,Tuple[float]]]

def main() -> None:
    '''Run functions assigned to CLI keywords or print help if no arguments passed.
    '''
    # get name of this module rather than __main__.py
    prog = __name__.rpartition('.')[-1]

    # setup incredibly basic logger
    logging.basicConfig(format='%(levelname)s: %(message)s')

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
    '''Gets assigned chemical shifts from BMRB entry, selects atoms, plots spectra
    and handes output.
    '''
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
    
    if not correlations:
        sys.exit(1)
    
    if args.csv is not None:
        write_csv(args.csv, correlations, selector.atoms, args.offset)

    _, ax = plt.subplots()
    # Plot correlations
    if len(selector.atoms) == 2:
        plotting.plot2D(ax, correlations, selector.atoms, color=args.color, nolabels=args.nolabels, showlegend=args.showlegend, offset=args.offset)

    elif len(selector.atoms) == 3:
        # Have to return Slices3D object, otherwise gets GC'd and scrolling doesn't work.
        three_d = plotting.plot3D(ax, correlations, selector.atoms, color=args.color, nolabels=args.nolabels, showlegend=args.showlegend, offset=args.offset,
        project=args.project-1, slices=args.slices)

    # Interactive matplotlib window opened if not saving
    if not args.output:
        plt.show()
    else:
        if len(selector.atoms) == 2:
            plt.savefig(args.output, dpi=300)
        else:
            # save each slice of 3D spectra as a separate file
            plt.close()
            atoms = selector.atoms[:args.project-1] + selector.atoms[args.project:]
            for i, data in enumerate(three_d.data):
                _, ax = plt.subplots()
                plotting.plot2D(ax, data, atoms, color=args.color, nolabels=args.nolabels, showlegend=args.showlegend, offset=args.offset)
                ax.annotate(f'{(sum(three_d.intervals[i]))/2:.1f} ppm', xy=(0.02,0.94), xycoords='axes fraction', horizontalalignment='left')

                name, extension = args.output.split('.')
                plt.savefig(f'{name}_slice{i+1}.{extension}', dpi=300)
                plt.close()
    plt.close()

def from_file(args: argparse.Namespace) -> None:
    '''Opens list of CSV files and plots overlay of all spectra.
    '''
    _, ax = plt.subplots()
    color_cycle = itertools.cycle(args.colors)
    all_correlations = []
    
     # need to keep track for showing CSPs
    csps = defaultdict(list)
    all_text = []
    fix_labels = False
    
    for in_file, color in zip(args.input, color_cycle):
        correlations = []
        with in_file as csvfile:
            reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
            atoms = tuple(next(reader))[1:] # get headings
            for row in reader:
                # sequence_number, residue_type, (shift1, shift2,...)
                correlations.append([int(row[0][3:]), row[0][:3], tuple(row[1:])])

        if all(len(correlation[-1]) == 2 for correlation in correlations):
            if args.showcsp and args.nolabels:
                # Use lables to match peaks for CSPs
                fix_labels = True
                args.nolabels = False
            _, text = plotting.plot2D(ax, correlations, atoms, color=color, nolabels=args.nolabels, showlegend=args.showlegend)
            for label in text:
                all_text.append(label)
                csps[label.get_text()].append(label.xy)
        else:
            all_correlations.extend(correlations)
    
    if all(len(correlation[-1]) == 3 for correlation in correlations):
        three_d = plotting.plot3D(ax, all_correlations, atoms, color=color, nolabels=args.nolabels, showlegend=args.showlegend,
        project=args.project-1, slices=args.slices)

    if args.showcsp:
        for points in csps.values():
            if len(points) > 1:
                xs, ys = zip(*points)
                plt.plot(xs, ys, 'k--', linewidth=1)
        if fix_labels:
            for label in all_text:
                label.remove()

    if args.csv is not None:
        write_csv(args.csv, all_correlations, atoms)

    # Interactive matplotlib window opened if not saving
    if not args.output:
        plt.show()
    else:
        if all(len(correlation[-1]) == 2 for correlation in correlations):
            plt.savefig(args.output, dpi=300)
        else:
            # save each slice of 3D spectra as a separate file
            plt.close()
            atoms = atoms[:args.project-1] + atoms[args.project:]
            for i, data in enumerate(three_d.data):
                _, ax = plt.subplots()
                plotting.plot2D(ax, data, atoms, color=color, nolabels=args.nolabels, showlegend=args.showlegend)
                ax.annotate(f'{(sum(three_d.intervals[i]))/2:.1f} ppm', xy=(0.02,0.94), xycoords='axes fraction', horizontalalignment='left')

                name, extension = args.output.split('.')
                plt.savefig(f'{name}_slice{i+1}.{extension}', dpi=300)
                plt.close()
    plt.close()

def search(args: argparse.Namespace) -> None:
    '''Wrapper for bmrb.instant_search which prints search results on the command line.
    '''
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
    '''Wrapper for bmrb.open_website.
    '''
    bmrb.open_website(args.entry)

def write_csv(outfile: TextIO, correlations: Correlations, atoms: Tuple[str], offset: int = 0) -> None:
    '''Helper function for outputting CSV files when --csv is passed.
    '''
    with outfile as csvfile:
        fieldnames = ['label'] + list(atoms)
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, quoting=csv.QUOTE_NONNUMERIC)
        writer.writeheader()
        for sequence_number, residue_type, chemical_shifts in correlations:
            row = {atom : chemical_shifts[i] for i, atom in enumerate(atoms)}
            row['label'] = f'{residue_type}{sequence_number + offset}'
            writer.writerow(row)

def choose_entity(names: list) -> int:
    '''Helper function for entries with more than one entity, prompts user with choice of which entity.
    '''
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
#!/usr/bin/env python

'''Produce simulated 2D spectra for protein entries from the BMRB using assigned chemical shifts.
This script provies functions for plotting commonly used correlations such as amide and methyl
as well as more complex custom correltations.

Peak assignments are pulled using the BMRB API (https://github.com/uwbmrb/BMRB-API) and plotted with matplotlib. 
Spectra produced by this script do not always reflect reality, as many experiments cannot create these correlations (yet!).
Many entries do not have 100% assignments so be sure to know what assignments are contained in the entry beforehand.
'''

import argparse
import itertools
import warnings
from collections import defaultdict
from operator import itemgetter

import requests
import matplotlib.pyplot as plt
import matplotlib.cm as cm

residue_map = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU',
               'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'M':'MET',
               'K':'LYS', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP',
               'Y':'TYR', 'V':'VAL'}
# Methyl atom names for MILVAT
methyl_atoms = {'ILE':[('HD11','CD1')], 'LEU':[('HD11','CD1'), ('HD21','CD2')],
                'VAL':[('HG11','CG1'), ('HG21','CG2')], 'MET':[('HE1','CE')],
                'ALA':[('HB1','CB')], 'THR':[('HG21','CG2')]}
simple_atoms = {'H', 'HA', 'N', 'C', 'CA', 'CB'}
# Palette of distinct colors for each residue
residue_colors = {'ALA':'maroon', 'ARG':'red', 'ASN':'pink', 'ASP':'brown',
                  'CYS':'orange', 'GLU':'coral', 'GLN':'olive', 'GLY':'magenta',
                  'HIS':'khaki', 'ILE':'purple', 'LEU':'green', 'MET':'navy',
                  'LYS':'blue', 'PHE':'lime', 'PRO':'lightgreen', 'SER':'aquamarine',
                  'THR':'cyan', 'TRP':'black', 'TYR':'grey', 'VAL':'yellow'}

class IndexTracker:
    # Modified from: https://matplotlib.org/stable/gallery/event_handling/image_slices_viewer.html
    def __init__(self, ax, data, intervals):
        self.ax = ax
        self.data = data
        self.intervals = intervals
        self.slices = len(self.data)
        self.ind = self.slices//2 # start in the center
        self.text = []
        self.handles = []
        try:
            for name, x, y in self.data[self.ind]:
                color = residue_colors[name[:3]]
                self.handles.append(self.ax.plot(x, y, 'o', c=color)[0])
                self.text.append(self.ax.annotate(name, (x,y)))
        except ValueError:
            # bin is empty, skip
            pass
        # Add ppm of center projected dimension
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        x_padding = 0.1 * (xmax - xmin)
        y_padding = 0.1 * (ymax - ymin)
        self.text.append(self.ax.text(xmin+0.1*x_padding, ymax-0.5*y_padding, f'{(sum(self.intervals[self.ind]))/2:.1f} ppm'))
        # self.update()

    def on_scroll(self, event):
        if event.button == 'up':
            # up goes down in ppm, backward because axes are backward
            if self.ind > 0:
                self.ind -= 1
                self.update()
        else:
            if self.ind < self.slices - 1:
                self.ind += 1
                self.update()

    def update(self):
        # Clear annotations
        try:
            for t in self.text:
                t.remove()
            self.text = []
            for h in self.handles:
                h.set_data(None,None)
            self.handles = []
            for name, x, y in self.data[self.ind]:
                color = residue_colors[name[:3]]
                self.handles.append(self.ax.plot(x, y, 'o', c=color)[0])
                self.text.append(self.ax.annotate(name, (x,y)))
        except ValueError:
            # bin is empty, clear plot
            for h in self.handles:
                h.set_data(None, None)
            for t in self.text:
                t.remove()
            self.text = []

        # Add ppm of center projected dimension
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        x_padding = 0.1 * (xmax - xmin)
        y_padding = 0.1 * (ymax - ymin)
        self.text.append(self.ax.text(xmin+0.1*x_padding, ymax-0.5*y_padding, f'{(sum(self.intervals[self.ind]))/2:.1f} ppm'))
        self.ax.figure.canvas.draw()

def get_bmrb_shifts(entry):
    # Use BMRB API to get entry number's entity names and assigned chemical shifts
    # Cannot requests loops and saveframes at the same time, so two requests
    names_url = f'http://api.bmrb.io/current/entry/{entry}?tag=Entity.name'
    shift_url = f'http://api.bmrb.io/current/entry/{entry}?loop=Atom_chem_shift'

    
    # Parses JSON returned by BMRB into a dictionary
    with requests.get(names_url) as name_tag:
        name_json = name_tag.json()
        # Handle mistyped entry codes or missing data
        if list(name_json.keys())[0] == 'error':
            raise ValueError(name_json['error'])
        entity_names = name_json[entry]['Entity.name']

    # Select which entity to plot if more than one
    entity_selection = 0
    if len(entity_names) > 1:
        while True:
            # Handle invalid inputs, break when valid
            try:
                print('Found more than one assigned chemical shift list, please select:')
                # Print list of names received from names_url request number matches Entity_ID tag
                for i, name in enumerate(entity_names, start=1):
                    print(f'({i}) {name}')
                entity_selection = int(input('> ')) - 1 # Entity_ID is 1-indexed, python is 0-indexed
                if 0 <= entity_selection < len(entity_names):
                    # Valid input
                    break
                else:
                    entity_selection = print(f'Invalid selection please choose a number <= {len(entity_names)}.')
            except ValueError:
                entity_selection = print(f'Invalid selection, please input a number <= {len(entity_names)}.')

    
    # Use BMRB API to get entry number's assigned chemical shifts
    with requests.get(shift_url) as shift_loop:
        bmrb_json = shift_loop.json()
        bmrb_dict = bmrb_json[entry]['Atom_chem_shift'][entity_selection]
        # Formats the data as a list where values are indexed by their NMRstar tag i.e. data_table[0]['Seq_ID']
        # This is a little bit more explicit than just using the column numbers
        data_table = [dict(zip(bmrb_dict['tags'], row)) for row in bmrb_dict['data']]
    return data_table

def get(args):
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
    plus_minus = None
    if args.custom is not None:
        args.custom = [a.upper() for a in args.custom]
        # This is the same as the methyl argument and this makes it easier to handle
        # if 'CMETHYL' in args.custom and 'HMETHYL' in args.custom:
        #     args.methyl = []
        #     args.custom = None
        for i, custom_atom in enumerate(args.custom):
            if '-' in custom_atom:
                if plus_minus is not None:
                    raise ValueError('only one atom can have +-')
                plus_minus = -1
                # Save which atom choice has the offset
                args.custom[i], mod_val = custom_atom.split('-')
                mod_val = int(mod_val)
                mod_index = i
            elif '+' in custom_atom:
                if plus_minus is not None:
                    raise ValueError('only one atom can have +-')
                plus_minus = 1
                # Save which atom choice has the offset
                args.custom[i], mod_val = custom_atom.split('+')
                mod_val = int(mod_val)
                mod_index = i
    if args.amide:
        # Make a dictionary similar to methyl_atoms
        selector = dict(zip(residues, [[('H', 'N')]]*len(residue_map.values())))
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

    data_table = get_bmrb_shifts(args.entry)

    # Parse BMRB assignments file storing relevant fields
    # Only use residue number, residue name, atom name and chemical shift
    
    assignments = defaultdict(dict)
    encountered_atoms = defaultdict(list)
    for row in data_table:
        # format for assignments dict {RES<num>: {atom1: shift, atom2: shift...}}
        assignments[row['Comp_ID']+row['Seq_ID']][row['Atom_ID']] = float(row['Val']) # Val is value, not valine
        # encountered_atoms is formatted with atom names indexed by residue type i.e. {'ALA':['H','N','C','CA'...]}
        if encountered_atoms.get(row['Comp_ID']) is None or row['Atom_ID'] not in encountered_atoms[row['Comp_ID']]:
            encountered_atoms[row['Comp_ID']].append(row['Atom_ID'])
    
    # Generate selectors for custom correlations from encountered atom names
    if args.custom is not None:
        selector_dicts = []
        for custom_atom in args.custom:
            # Hmethyl and Cmethyl already defined by modifying methyl_atoms
            if custom_atom == 'HMETHYL':
                h_selector = {}
                for res, pairs in methyl_atoms.items():
                    h_selector[res] = [p[0] for p in pairs]
                selector_dicts.append(h_selector)
            elif custom_atom == 'CMETHYL':
                c_selector = {}
                for res, pairs in methyl_atoms.items():
                    c_selector[res] = [p[1] for p in pairs]
                selector_dicts.append(c_selector)
            # If the atom is not a simple case match against the first letters of the atom IDs
            elif custom_atom not in simple_atoms:
                custom_selector = {}
                for res, atoms in encountered_atoms.items():
                    custom_selector[res] = [a for a in atoms if a.startswith(custom_atom)]
                selector_dicts.append(custom_selector)
            # If custom_atom is the simple case, generate dict similar to how amide is done
            else:
                # residue_map.values() contains all residues, does not apply filter yet 
                simple_selector = dict(zip(residue_map.values(), [[(custom_atom)]]*len(residue_map.values())))
                selector_dicts.append(simple_selector)
        
        selector = {}
        if plus_minus is not None:
            # will need selectors to be based on residue pairs {res1_res2: [(atom1, atom2), ...]}
            key_groups = itertools.product(*(sd.keys() for sd in selector_dicts))
            for keys in key_groups:
                groups = (selector_dicts[i][res] for i,res in enumerate(keys))
                selector['_'.join(keys)] = list(itertools.product(*groups))
        else:
        # Does not matter which is longer, will either be filtered here or when plotting
            for res in residue_map.values():
                try:
                    groups = (sd[res] for sd in selector_dicts)
                    selector[res] = list(itertools.product(*groups))
                except KeyError:
                    # Key was not present in all dicts, ignore
                    continue

    # Data structure containing selected atoms' peaks
    resonances = defaultdict(dict)

    # For an i+- correlation
    if plus_minus is not None:
        # Split into two iterators high is larger residue number, low is the smaller one
        low, high = itertools.tee(assignments.items(), 2)
        # Shift the high iterator forward by mod_val
        high = itertools.islice(high, mod_val, None)

        # keep the i residue in the second position, swapped for + and -
        if plus_minus < 0:
            split_iterator = itertools.zip_longest(low, high, fillvalue=None)
        elif plus_minus > 0:
            split_iterator = itertools.zip_longest(high, low, fillvalue=None)

        # Either first or second index
        ires_index = 0 if mod_index == 1 else 1
        for mod_assign, ires_assign in split_iterator:
            try:
                # first index in iterator is residue name/num second is atom dict
                mod_atom_dict = mod_assign[1]
                ires_atom_dict = ires_assign[1]
                ires = ires_assign[0][:3] # residue type of i-residue
                mod_res = mod_assign[0][:3] # residue type of +- residue

                # Check for gaps in the sequence, plus_minus is +-1
                if int(ires_assign[0][3:]) + plus_minus * mod_val == int(mod_assign[0][3:]):
                    if ires_index == 0:
                        pairs = selector.get(f'{ires}_{mod_res}')
                    else:
                        pairs = selector.get(f'{mod_res}_{ires}')
                    res_tuple = [None] * len(args.custom)
                    # Apply reisude filter
                    if ires in residues and pairs is not None:
                        for p in pairs:
                            # Insert data into proper location of tuple based which atom has +-
                            res_tuple[ires_index] = ires_atom_dict.get(p[ires_index])
                            res_tuple[mod_index] = mod_atom_dict.get(p[mod_index])
                            if ires_index == 0:
                                if args.label == 1:
                                    resonances[ires_assign[0]]['_'.join(p)] = tuple(res_tuple)
                                else:
                                    resonances[mod_assign[0]]['_'.join(p)] = tuple(res_tuple)
                            else:
                                if args.label == 1:
                                    resonances[mod_assign[0]]['_'.join(p)] = tuple(res_tuple)
                                else:
                                    resonances[ires_assign[0]]['_'.join(p)] = tuple(res_tuple)
            except TypeError:
                # Reached None the end of the curr iterator filled with None by zip_longest
                continue
    else:
    # For a typical spectrum
        for name, atom_dict in assignments.items():
            res = name[:3]
            groups = selector.get(res)
            # Apply residue filter
            if res in residues and groups is not None:
                for atoms in groups:
                    # Uses get to fill None if a chemical shift is not present
                    resonances[name]['_'.join(atoms)] = tuple(atom_dict.get(a) for a in atoms)

    seq_nums = [int(idx[3:]) for idx in resonances.keys()]
    seq_min = min(seq_nums)
    seq_max = max(seq_nums)
    if args.segment is not None:
        start, stop = args.segment
        if start < seq_min:
            warnings.warn(f'Start from segment parameter: {start} too low, defaulting to sequence start: residue {seq_min}')
            start = seq_min 
        if stop > seq_max:
            warnings.warn(f'Stop from segment parameter: {stop} too high, defaulting to sequence end: residue {seq_max}')
            stop = seq_max
        if start > seq_max:
            warnings.warn(f'Start from segment parameter: {stop} too high, ignoring')
            start = seq_min

    legend = {}
    
    # 3D plotting
    if args.custom is not None and len(args.custom) == 3:
        plot_list = []
        # Reformat into a list of lists similar to CSV file [[name, peak1, peak2, peak3],...]
        for name, pairs in resonances.items():
            idx = int(name[3:])
            # Only keep specfied segment if segment argument is passed
            if args.segment is not None and not start <= idx <= stop:
                continue
            for atoms, peak in pairs.items():
                if None not in peak:
                    plot_list.append([name] + list(peak))
        
        # Get mins and maxes so all plots have the same xy coords
        maxes = [max(dim) for dim in list(zip(*plot_list))[1:]]
        mins = [min(dim) for dim in list(zip(*plot_list))[1:]]
        
        # Find intervals for bins
        proj_max = maxes.pop(args.project-1)
        proj_min = mins.pop(args.project-1)
        bin_width = (proj_max - proj_min) / args.slices
        cutoffs = [proj_min + i * bin_width for i in range(args.slices)]
        # Add 1 to not have to deal with < vs <= for last bin
        intervals = list(zip(cutoffs, cutoffs[1:] + [proj_max + 1]))

        fig, ax = plt.subplots()
        # Set axes equal for all plots
        x_padding = 0.1 * (maxes[0] - mins[0])
        y_padding = 0.1 * (maxes[1] - mins[1])
        ax.set_xlim((mins[0]-x_padding, maxes[0]+x_padding))
        ax.set_ylim((mins[1]-y_padding, maxes[1]+y_padding))
        
        # Invert axes
        ax.invert_xaxis()
        ax.invert_yaxis()

        sliced_data = []
        for low, high in intervals:
            sliced = []
            for point in sorted(plot_list, key=itemgetter(args.project)):
                res = point[0][:3]
                idx = int(point[0][3:])
                if low <= point[args.project] < high:
                    # Remove projected index from list of indices used for plotting
                    plot_point = point[:]
                    plot_point.pop(args.project)
                    sliced.append([f'{res}{idx+args.offset}'] + plot_point[1:])
                    # Plot 2D projections
                    # handle = ax.scatter(*plot_point[1:], c=residue_colors[res])
                    # legend[res] = handle
                    # Show the average of the bin interval in the upper center
                    # if not args.nolabels:
                    #     ax.annotate(f'{res}{idx+args.offset}', plot_point[1:])
            sliced_data.append(sliced)

        tracker = IndexTracker(ax, sliced_data, intervals)
        fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
        plt.show()
    
    else:
        # Plottting
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
            plt.xlabel(f'{args.custom[0]} - {nucleus1} (ppm)')
            plt.ylabel(f'{args.custom[1]} - {nucleus2} (ppm)')

        # Create output file
        if args.csv:
            output = open(args.csv, 'w')

        legend = {}
        for name, pairs in resonances.items():
            res = name[:3]
            idx = int(name[3:])
            # Only plot specfied segment if segment argument is passed
            if args.segment is not None and not start <= idx <= stop:
                continue
            for atoms, peak in pairs.items():
                if None not in peak:
                    # Plot peak centers colored by residue
                    handle = plt.scatter(*peak, c=residue_colors[res])
                    legend[res] = handle
                    # Label each peak, add atom ID to prochiral atoms for Leu/Val
                    if (args.methyl is not None and res in {'LEU', 'VAL'}):
                        # Adds the atom name (i.e. CD1 to the end of each label)
                        peak_label = f"{res}{idx+args.offset}{atoms.split('_')[1]}"
                        if not args.nolabels:
                            plt.annotate(peak_label, peak)
                        if args.csv:
                            output.write(f'{peak_label},' + ','.join(str(p) for p in peak) + '\n')
                    else:
                        # Label is residue type and sequence index
                        if not args.nolabels:
                            plt.annotate(f'{res}{idx+args.offset}', peak)
                        if args.csv:
                            output.write(f'{res}{idx+args.offset},' + ','.join(str(p) for p in peak) + '\n')

        # Show the legend when no labes are shown, or showlegend argument is passed
        if args.nolabels or args.showlegend:
            # Residues are added as they are encountered, not in sorted order
            sort_legend = {k: v for k, v in sorted(legend.items())}
            plt.legend(sort_legend.values(), sort_legend.keys())

        # Invert axes
        plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()

        if args.csv:
            output.close()

        # Interactive matplotlib window opened if not saving
        if not args.output:
            plt.show()
        else:
            plt.savefig(args.output, dpi=300)
        plt.close()

def plot(args):
    # Create output file
    if args.csv:
        output = open(args.csv, 'w')

    legend = {}
    for csv in args.input:
        with open(csv, 'r') as f:
            shifts = [line.strip().split(',') for line in f]
        for name, x, y in shifts:
            res = name[:3]
            color = residue_colors[res]
            handle = plt.scatter(float(x), float(y), c=color)
            legend[res] = handle
            
            if not args.nolabels:
                plt.annotate(name, (float(x),float(y)))
            if args.csv:
                output.write(f'{name},{x},{y}\n')
        
        # Show the legend when no labes are shown, or showlegend argument is passed
        if args.nolabels or args.showlegend:
            # Residues are added as they are encountered, not in sorted order
            sort_legend = {k: v for k, v in sorted(legend.items())}
            plt.legend(sort_legend.values(), sort_legend.keys())

    # Invert axes
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    
    if args.csv:
        output.close()

    # Interactive matplotlib window opened if not saving
    if not args.output:
        plt.show()
    else:
        plt.savefig(args.output, dpi=300)
    plt.close()

def search(args):
    search_url = f'http://api.bmrb.io/current/instant?term={args.query}&database=macromolecules'
    with requests.get(search_url) as search_results:
        for results in search_results.json():
            data = set(r.get('type') for r in results['data_types'])
            if 'assigned_chemical_shifts' in data:
                print(f"{results['value']}\t{results['label']}")

def main():
    # CLI flags with descriptions
    parser = argparse.ArgumentParser(description='Tools for plotting simulated NMR spectra from assigned chemical shifts in the BMRB.')
    subparsers = parser.add_subparsers()

    # Subcommand `get` is the primary function of the script
    get_parser = subparsers.add_parser('get', description='Pull assigned chemical shifts from BMRB entry and plot peak positions on a simulated 2D spectrum.')
    get_parser.add_argument('entry', help='BMRB accession number')
    get_parser.add_argument('-r', '--residues', help='Residues filter: a single string (i.e. ILV) of 1-letter amino acid codes. Will only plot these residues.')
    get_parser.add_argument('--offset', default=0, type=int, help='Index offset. BMRB indices starts at 1 which may not reflect the true residue number')
    get_parser.add_argument('--segment', nargs=2, type=int, help='Indices representing the start and stop of the sequence you want to show. Does not account for offset Format: start stop')
    get_parser.add_argument('--csv', help='Name for output file of shifts in CSV format for later plotting.')
    get_parser.add_argument('--label', default=1, type=int, help='Choose which residue is used for the label in custom correlations: 1, 2 or 3. Default is 1.')
    get_parser.add_argument('--project', type=int, help='Choose which dimension to project on for 3D spectra, 1, 2 or 3.')
    get_parser.add_argument('--slices', default=16, type=int, help='Number of slices/points in 3rd dimension in 3D spectra.')
    get_parser.add_argument('--nolabels', action='store_true', help='Plot spectrum without peak labels.')
    get_parser.add_argument('--showlegend', action='store_true', help='Show a legend on the spectrum with residue colors.')
    get_parser.add_argument('-o', '--output', help='Save spectrum to file and do not preview plot with matplotlib.')

    # Flags to set which correlation will be depicted
    correlation = get_parser.add_mutually_exclusive_group(required=True)
    correlation.add_argument('--amide', action='store_true', help='Plot amide shifts.')
    correlation.add_argument('--methyl', nargs='*', choices=['proR', 'proS'], help='Plot methyl shifts, allows use of proR and proS flags for specific labeling schemes. By default plots both.')
    correlation.add_argument('--custom', nargs='*', help='Plot correlation between any two atoms (i.e. N CO).')
    get_parser.set_defaults(func=get)

    # Allows spectra saved to CSV files to be replotted or superimposed
    plot_parser = subparsers.add_parser('plot', description='Plot one or more spectra from CSV output files saved by `get`. Can also output combined CSVs for multiple spectra.')
    plot_parser.add_argument('input', nargs='*', help='Name of input CSV files to overlay.')
    plot_parser.add_argument('--csv', help='Name for output file of shifts in CSV format for later plotting.')
    plot_parser.add_argument('--nolabels', action='store_true', help='Plot spectrum without peak labels.')
    plot_parser.add_argument('--showlegend', action='store_true', help='Show a legend on the spectrum with residue colors.')
    plot_parser.add_argument('-o', '--output', help='Save spectrum to file and do not preview plot with matplotlib.')
    plot_parser.set_defaults(func=plot)

    # Search BMRB for an entry ID using text or a PDB ID
    search_parser = subparsers.add_parser('search', description='Search the BMRB for an entry number for use with `get`. Can use simple search terms including a PDB code')
    search_parser.add_argument('query', nargs='?', help='Search terms for BMRB database.')
    search_parser.set_defaults(func=search)

    # Parse arguments and run assigned function
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
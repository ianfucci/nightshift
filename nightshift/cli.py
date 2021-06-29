import argparse

import matplotlib.pyplot as plt

import nightshift.nightshift

def build_parser(prog: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description='''A python command line program for plotting simulated 2D and 3D
            NMR spectra from assigned chemical shifts from the BMRB''',
        epilog='''Use -h or --help with each command to get help (i.e. %(prog)s get -h)''',
        )
    subparsers = parser.add_subparsers()

    # Subcommand `get` is the primary function of the script
    get_parser = subparsers.add_parser('get', 
        description='''Get assigned chemical shifts from BMRB entry and plot peak 
            positions on a simulated spectrum''',
        help='''Get assigned chemical shifts from BMRB entry and plot peak 
            positions on a simulated spectrum''',
        )
    get_parser.add_argument('entry',
        help='''BMRB entry number. If not known, try `%(prog)s search`''',
        )

    # Set which correlation will be depicted
    correlation = get_parser.add_mutually_exclusive_group(required=True)
    correlation.add_argument('--amide',
        action='store_true',
        help='''Plot simulated 1H-15N HSQC spectrum''',
        )
    correlation.add_argument('--methyl',
        action='store_true',
        help='''Plot simulated 1H-13C HMQC spectrum, allows use of --proR and --proS
            for specific labeling schemes''',
        )
    correlation.add_argument('--custom',
        nargs='*',
        metavar=('ATOM1', 'ATOM2'),
        help='''Plot correlation between any two or three atoms (i.e. N C).
            Including +<number> or -<number> after an atom name allows for inter-residue
            correlations (i.e. H N C-1). Also supports HMETHYL and CMETHYL and proR/proS.''',
        )
    
    get_parser.add_argument('-r', '--residues',
        help='''Residue filter is  single string of 1-letter amino acid codes
            (i.e. ILV). Will select only these residues''',
        )
    get_parser.add_argument('--offset',
        default=0,
        type=int,
        help='''Add a value to BMRB index numbers, most entries start with the first
            residue as 1, this may not be true (default: %(default)s)'''
        )
    get_parser.add_argument('--segment',
        nargs=2,
        type=int,
        metavar=('START', 'STOP'),
        help='''Indices for the residue numbers and stop of the sequence you want to show,
             uses indices from BMRB file and does not account for offset''',
        )
    get_parser.add_argument('--csv',
        metavar='FILENAME', 
        help='''Name for output file of shifts in CSV format for later plotting'''
        )
    get_parser.add_argument('--label', 
        type=int,
        choices=[1,2,3],
        help='''Choose which residue is used for the label in custom correlation
            (choices: %(choices)s) (default: i residue)'''
        )
    get_parser.add_argument('--project', 
        type=int, 
        choices=[1,2,3],
        default=3, 
        help='''Choose which dimension to project on for 3D spectra (choices: %(choices)s)
            (default: %(default)s)'''
        )
    get_parser.add_argument('--slices', 
        default=16,
        type=int,
        help='''Number of slices/points in projected dimension of a 3D spectrum
            (default: %(default)s)''',
        )
    get_parser.add_argument('--nolabels', 
        action='store_true',
        help='''Plot spectrum without peak labels''',
        )
    get_parser.add_argument('--showlegend', 
        action='store_true', 
        help='''Show a legend on the spectrum with residue colors''',
        )
    
    mpl_filetypes = plt.gcf().canvas.get_supported_filetypes()
    file_extensions = ' '.join(f'.{extension}' for extension in mpl_filetypes.keys())
    get_parser.add_argument('-o', '--output', 
        metavar='FILENAME',
        help=f'''Save spectrum to file and do not plot in interactive mode with matplotlib.
            Can have the following extensions: {file_extensions} '''
        )

    # Show sidechains for --amide or proR/proS for --methyl and --custom HMETHYL/CMETHYL
    atom_filters = get_parser.add_mutually_exclusive_group()
    atom_filters.add_argument('--sidechains', 
        action='store_true', 
        help='''Choose to plot side chains of Asn and Gln, use with --amide''',
        )
    atom_filters.add_argument('--proR', 
        action='store_true', 
        help='''Choose to plot only proR peaks for Leu and Val, use with --methyl
            and --custom''',
            )
    atom_filters.add_argument('--proS', 
        action='store_true', 
        help='''Choose to plot only proS peaks for Leu and Val, use with --methyl
            and --custom''',
            )
    
    # Allows spectra saved to CSV files to be replotted or superimposed
    open_parser = subparsers.add_parser('open',
        description='''Plot one or more spectra from CSV output files saved by `get`.
            Allows spectra saved to CSV files to be replotted or superimposed.''',
        help='''Plot one or more spectra from CSV output files saved by `get`.
            Allows spectra saved to CSV files to be replotted or superimposed.''',
        )
    open_parser.add_argument('input',
        nargs='+',
        metavar='FILENAME(S)',
        help='''Name(s) of input CSV file(s) to plot.''',
        )
    open_parser.add_argument('--csv',
        metavar='FILENAME', 
        help='''Name for output file of shifts in CSV format for later plotting'''
        )
    open_parser.add_argument('--label', 
        type=int,
        choices=[1,2,3],
        help='''Choose which residue is used for the label in custom correlation
            (choices: %(choices)s) (default: i residue)'''
        )
    open_parser.add_argument('--project', 
        type=int, 
        choices=[1,2,3],
        default=3, 
        help='''Choose which dimension to project on for 3D spectra (choices: %(choices)s)
            (default: %(default)s)'''
        )
    open_parser.add_argument('--slices', 
        default=16,
        type=int,
        help='''Number of slices/points in projected dimension of a 3D spectrum
            (default: %(default)s)''',
        )
    open_parser.add_argument('--nolabels', 
        action='store_true',
        help='''Plot spectrum without peak labels''',
        )
    open_parser.add_argument('--showlegend', 
        action='store_true', 
        help='''Show a legend on the spectrum with residue colors''',
        )
    open_parser.add_argument('-o', '--output', 
        metavar='FILENAME',
        help=f'''Save spectrum to file and do not plot in interactive mode with matplotlib.
            Can have the following extensions: {file_extensions} '''
        )

    # Search BMRB for an entry ID using text or a PDB ID
    search_parser = subparsers.add_parser('search',
        description='''Search the BMRB for an entry number for use with `get`.
        Can use simple search terms including a PDB code.''',
        help='''Search the BMRB for an entry number for use with `get`.
        Can use simple search terms including a PDB code.''',
        )
    search_parser.add_argument('terms', nargs='*', help='Search terms for BMRB database.')

    # Open BMRB webpage for an entry number
    website_parser = subparsers.add_parser('website', 
        description='''Open BMRB webpage for an entry number''',
        help='''Open BMRB webpage for an entry number''',
        )
    website_parser.add_argument('entry',
        help='''BMRB entry number. If not known, try `%(prog)s search`''',
        )
    
    # Assign functions to pass arguments to
    get_parser.set_defaults(func=nightshift.nightshift.get)
    open_parser.set_defaults(func=nightshift.nightshift.from_file)
    search_parser.set_defaults(func=nightshift.nightshift.search)
    website_parser.set_defaults(func=nightshift.nightshift.website)

    return parser
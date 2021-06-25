# nightshift
A python library and command line program for plotting simulated 2D and 3D NMR spectra from assigned chemical shifts from the [BMRB](https://bmrb.io)

Peak assignments are pulled using the BMRB API (https://github.com/uwbmrb/BMRB-API) and plotted with matplotlib. Spectra produced by this script do not always reflect reality, as many entries do not have 100% assignments. Be sure to know what assignments are contained in the entry beforehand.

## Installation
I'll assume the reader may not necessarily familiar with python and installing dependencies.
These scripts require python 3.7 or later, matplotlib and requests. I recommend installing in a virtual environment (so that you don't mess up anything else on your system). These instructions are adapted from https://docs.python.org/3/tutorial/venv.html

Navigate to the directory containing these scripts and run:
`python3 -m venv plot_bmrb`

This will create a new directory called plot_bmrb.
On Windows activate the virtual environment by running:

`plot_bmrb/Scripts/activate.bat`

or on Linux/macOS:

`source plot_bmrb/bin/activate`

_Note_: this virtual environment should be re-activated whenever you want to run these scripts.

Your shell should reflect that you have activated the virtual environment. You can check by running `which python` which should return that python points within this directory.

To install the dependencies run:
`pip install -r requirements.txt`

## Usage
You will need to know your protein of interest's BMRB entry number.
We'll use `4493` [Solution structure of the designed hydrophobic core mutant of ubiquitin, 1D7](http://www.bmrb.wisc.edu/data_library/summary/index.php?bmrbId=4493) for our examples.

### Basic spectra

Two major use cases are getting an idea what a amide or methyl spectrum of your protein would look like

To simulate 1H-15N HSQC spetra can be plotted using the `--amide` flag

`python bmrb_shifts.py get 4493 --amide`

Simulate 1H-13C HMQC spetra can be plotted using the `--methyl` flag, additionally
providing the optional `proS` or `proR` will filter LV atoms by prochirality

`python bmrb_shifts.py get 4493 --methyl proS`

All plots can be filtered based on residue type by passing the `--residues` or `-r` flag and providing one-letter amino acid codes. For instance an ILV methyl labeled spectrum can be plotted using

`python bmrb_shifts.py get 4493 --methyl -r ILV`

### Custom correlations

For arbitrary correlations use the `--custom` flag followed by two atom names. Consider yourself warned that labeling schemes and/or experiments to produce these correlations may not (currently) exist. Atoms for custom correlations are specified using standard PDB atom names: H for amide proton, N for amide nitrogen, C for carbonyl carbon, CA for alpha carbon, HA for alpha proton and so on. FOr particular residues two or more atoms may exist at a position (i.e. CG for Val could be CG1 or CG2). To specify both CG1 and CG2 for Val pass CG

`python bmrb_shifts.py get 4493 --custom CG CA -r V`

or specify the full atom name to only get those atoms

`python bmrb_shifts.py get 4493 --custom CG1 CA -r V`

Two special atom names also exist for custom correlations: `Hmethyl` and `Cmethyl`. Which correspond to these atoms in MILVAT residues and are the same atoms selected by using the `--methyl` flag.

| Residue | Hmethyl    | Cmethyl |
|---------|------------|---------|
| Met     | HE1        | CE      |
| Ile     | HD11       | CD1     |
| Leu     | HD11, HD21 | CD1,CD2 |
| Val     | HG11, HG21 | CG1,CG2 |
| Ala     | HB1        | CB      |
| Thr     | CG2        | HG21    |

This allows for correlations of methyl groups to any other atom to be plotted. For instance Cmethyl to CA

`python bmrb_shifts.py get 4493 --custom Cmethyl CA -r ILV`

Adding '-1' to the end of a custom atom name allows correlation to the _i_-1 residue. For instance correlation of the CO of the _i_-1 residue to the amide N of the _i_ residue

`python bmrb_shifts.py get 4493 --custom C-1 N`

### Other options

By default plots are generated in matplotlib and are interactive. To save directly to an image file use the `--output` or `-o` flag and provide a file name and extension (.eps, .pdf, .pgf, .png, .ps, .raw, .rgba, .svg, and .svgz are all acceptable).

Formatting options include:
- `--showlegend` to add a legend
- `--nolabels` to remove the residue/atom name and numbers from the plot, also shows the legend
- `--offset` to add a constant to the indices used by BMRB (to reflect the numbering you are used to)

A csv file containing the label and chemical shifts of both atoms can be saved using the `--csv` flag and providing a file name

`python bmrb_shifts.py get 4493 --methyl -r ILV --csv output.csv`

this can be opened in other software to generate plots with different formatting.
Also the auxiliary script `plot_outputs.py` can be used to overlay multiple spectra (i.e. different domains of the same protein or a protein complex).

This example is nonsense, but illustrates how it could be done. 
First, generate two output files

`python bmrb_shifts.py get 4493 --amide --csv output1.csv`

`python bmrb_shifts.py 3433 --amide --csv output2.csv`

Then plot their overlay (_Note_: plot_output.py is used here)

`python plot_output.py output1.csv output2.csv`

This script also accepts the `--showlegend` and `--nolabels` flags.

## Interesting examples
- ILV methyl spectrum:

  `python bmrb_shifts.py get 4493 --methyl -r ILV`

- ILV methyl spectrum with proR LV:

  `python bmrb_shifts.py get 4493 --methyl proR -r ILV`

- Amide spectrum of only lysines and arginines:

  `python bmrb_shifts.py get 4493 --amide -r KR`

- 2D HMBC-HMQC (intra-residue methyl-methyl correlations):

  `python bmrb_shifts.py get 4493 --custom Cmethyl Cmethyl -r LV`

- 2D NCO:

  `python bmrb_shifts.py get 4493 --custom C-1 N`

- Arg/Lys side chain carbon correlations (a la [Pritchard and Hansen, 2019](https://doi.org/10.1038/s41467-019-09743-4))
  
  `python bmrb_shifts.py get 4493 --custom CG CD -r R --csv 4493_arg.csv`
  
  `python bmrb_shifts.py get 4493 --custom CD CE -r K --csv 4493_lys.csv`
  
  `python bmrb_shifts.py plot 4493_arg.csv 4493_lys.csv --showlegend`
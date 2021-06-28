# nightshift
A python library and command line program for plotting simulated 2D and 3D NMR spectra from assigned chemical shifts from the [BMRB](https://bmrb.io)

Peak assignments are pulled using the BMRB API (https://github.com/uwbmrb/BMRB-API) and plotted with matplotlib. Spectra produced by this script do not always reflect reality, as many entries do not have 100% assignments. Be sure to know what assignments are contained in the entry beforehand.

## Installation
I recommend installing in a [virtualenv](https://docs.python.org/3/tutorial/venv.html) to avoid any conflicts with your python installation.
nightshift is available on the [PyPI](https://pypi.org/project/nightshift/) and can be installed with `pip`:

`pip install nightshift`

To install a local version downloaded from GitHub, use:

`pip install /path/to/nightshift`

## Usage
You will need to know your protein of interest's BMRB entry number. You can use nightshift's search command to find entry numbers:

`nightshift search ubiquitin`

We'll use `4493` [Solution structure of the designed hydrophobic core mutant of ubiquitin, 1D7](http://www.bmrb.wisc.edu/data_library/summary/index.php?bmrbId=4493) for our examples.

### Basic spectra

Two major use cases are getting an idea what a amide or methyl spectrum of your protein would look like

To simulate 1H-15N HSQC spetra can be plotted using the `--amide` flag:

`nightshift get 4493 --amide`

To show Asn and Gln side chain amides on the spectrum pass the optional `--sidechains` flag:

`nightshift get 4493 --amide --sidechains`

A simulated 1H-13C HMQC spetra can be plotted using the `--methyl` flag, additionally
providing the optional `proS` or `proR` will filter LV atoms by prochirality:

`nightshift get 4493 --methyl --proS`

All plots can be filtered based on residue type by passing the `--residues` or `-r` flag and providing one-letter amino acid codes. For instance an ILV methyl labeled spectrum can be plotted using:

`nightshift get 4493 --methyl -r ILV`

### Custom correlations

For arbitrary correlations use the `--custom` flag followed by two atom names. Consider yourself warned that labeling schemes and/or experiments to produce these correlations may not (currently) exist. Atoms for custom correlations are specified using standard PDB atom names: H for amide proton, N for amide nitrogen, C for carbonyl carbon, CA for alpha carbon, HA for alpha proton and so on. For particular residues two or more atoms may exist at a position (i.e. CG for Val could be CG1 or CG2). To specify both CG1 and CG2 for Val pass CG:

`nightshift get 4493 --custom CG CA -r V`

or specify the full atom name to only get those atoms:

`nightshift get 4493 --custom CG1 CA -r V`

Two special atom names also exist for custom correlations: `Hmethyl` and `Cmethyl`. Which correspond to these atoms in MILVAT residues and are the same atoms selected by using the `--methyl` flag. The flags `--proS` and `--proR` can be used with `Hmethyl` and `Cmethyl`

| Residue | Hmethyl    | Cmethyl |
|---------|------------|---------|
| Met     | HE1        | CE      |
| Ile     | HD11       | CD1     |
| Leu     | HD11, HD21 | CD1,CD2 |
| Val     | HG11, HG21 | CG1,CG2 |
| Ala     | HB1        | CB      |
| Thr     | CG2        | HG21    |

This allows for correlations of methyl groups to any other atom to be plotted. For instance Cmethyl to CA:

`nightshift get 4493 --custom Cmethyl CA -r ILV`

Adding '-' or '+' and any number to the end of a custom atom name allows correlation to the _i_+/- _num_ residue. For instance correlation of the CO of the _i_-1 residue to the amide N of the _i_ residue:

`nightshift get 4493 --custom C-1 N`

### 3D correlations

The `--custom` option also allows for 3D correlations to be plotted. The matplotlib window will show 2D slices of the simulated spectrum. Scrolling the mouse wheel will switch between each slice:

`nightshift get 4493 --custom H N CA`

By default 16 slices are generated, this can be altered with the `--slices` option (i.e. `slices 32` or `slices 1` for a 2D projection). The `--project` parameter can take a value of 1, 2 or 3 which chooses which dimension to project on. For a 3D HNCA (though this can currently only be for i or i-1):

`nightshift get 4493 --custom H N CA --project 2 --slices 32`

Plus and minus can also be used on 3D correlations:

`nightshift get 4493 --custom HA N+1 CA+2 --project 1`

### Other options

By default plots are generated in matplotlib and are interactive. To save directly to an image file use the `--output` or `-o` flag and provide a file name and extension (.eps, .pdf, .pgf, .png, .ps, .raw, .rgba, .svg, and .svgz are all acceptable).

Formatting options include:
- `--showlegend` to add a legend
- `--nolabels` to remove the residue/atom name and numbers from the plot, also shows the legend
- `--offset` to add a constant to the indices used by BMRB (to reflect the numbering you are used to)

A csv file containing the label and chemical shifts of both atoms can be saved using the `--csv` flag and providing a file name:

`nightshift get 4493 --methyl -r ILV --csv output.csv`

this can be opened in other software to generate plots with different formatting.
Also the auxiliary script `plot_outputs.py` can be used to overlay multiple spectra (i.e. different domains of the same protein or a protein complex).

This example is nonsense, but illustrates how it could be done. 
First, generate two output files:

`nightshift get 4493 --amide --csv output1.csv`

`nightshift 3433 --amide --csv output2.csv`

Then plot their overlay:

`nightshift open output1.csv output2.csv`

This script also accepts the `--showlegend` and `--nolabels` flags.

## Interesting examples
- ILV methyl spectrum:

  `nightshift get 4493 --methyl -r ILV`

- ILV methyl spectrum with proR LV:

  `nightshift get 4493 --methyl proR -r ILV`

- Amide spectrum of only lysines and arginines:

  `nightshift get 4493 --amide -r KR`

- 2D HMBC-HMQC (intra-residue methyl-methyl correlations):

  `nightshift get 4493 --custom Cmethyl Cmethyl -r LV`

- 2D NCO:

  `nightshift get 4493 --custom C-1 N`

- Arg/Lys side chain carbon correlations (a la [Pritchard and Hansen, 2019](https://doi.org/10.1038/s41467-019-09743-4))
  
  `nightshift get 4493 --custom CG CD -r R --csv 4493_arg.csv`
  
  `nightshift get 4493 --custom CD CE -r K --csv 4493_lys.csv`
  
  `nightshift open 4493_arg.csv 4493_lys.csv --showlegend`
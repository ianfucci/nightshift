#!/usr/bin/env python

'''Auxiliary script for plotting output files from bmrb_shifts with the same format.
Used for taking multiple output files and plotting on the same spectrum.
'''

import sys
import matplotlib.pyplot as plt

# Palette of distinct colors for each residue from bmrb_shifts
residue_colors = {'ALA':'maroon', 'ARG':'red', 'ASN':'pink', 'ASP':'brown',
                  'CYS':'orange', 'GLU':'coral', 'GLN':'olive', 'GLY':'magenta',
                  'HIS':'khaki', 'ILE':'purple', 'LEU':'green', 'MET':'navy',
                  'LYS':'blue', 'PHE':'lime', 'PRO':'lightgreen', 'SER':'aquamarine',
                  'THR':'cyan', 'TRP':'black', 'TYR':'grey', 'VAL':'yellow'}

legend = {}
for csv in sys.argv[1:]:
    try:
        with open(csv, 'r') as f:
            shifts = [line.strip().split(',') for line in f]
        for name, x, y in shifts:
            res = name[:3]
            color = residue_colors[res]
            handle = plt.scatter(float(x), float(y), c=color)
            legend[res] = handle
            if '--nolabels' not in sys.argv:
                plt.annotate(name, (float(x),float(y)))
        if '--nolabels' in sys.argv or '--showlegend' in sys.argv:
            sort_legend = {k: v for k, v in sorted(legend.items())}
            plt.legend(sort_legend.values(), sort_legend.keys())
    except FileNotFoundError:
        if csv == '--nolabels' or csv == '--showlegend':
            pass
        else:
            print(f'"{csv}" not found, skipping.')
# Invert axes
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

plt.show()
plt.close()
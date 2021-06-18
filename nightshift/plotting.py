from typing import Dict, List, NamedTuple, Tuple
from operator import itemgetter

import numpy as np
import matplotlib.pyplot as plt


# Palette of colors for each residue
RESIDUE_COLORS: Dict[str,str] = {
    'ALA':'maroon', 'ARG':'red', 'ASN':'pink', 'ASP':'brown',
    'CYS':'orange', 'GLU':'coral', 'GLN':'olive', 'GLY':'magenta',
    'HIS':'khaki', 'ILE':'purple', 'LEU':'green', 'MET':'navy',
    'LYS':'blue', 'PHE':'lime', 'PRO':'lightgreen', 'SER':'aquamarine',
    'THR':'cyan', 'TRP':'black', 'TYR':'grey', 'VAL':'yellow',
    }

AXIS_LABELS: Dict[str,str] = {'H': r'$^1$H (ppm)', 'N': r'$^{15}$N (ppm)', 'C': r'$^{13}$C (ppm)'}

def plot2D(records: List[NamedTuple], *, nolabels: bool, showlegend: bool, offset: int, segment: Tuple[int]) -> None:
    # should really take care of segment and offset elsewhere... selector probably
    legend = {}
    fig, ax = plt.subplots()
    for record in records:
        index = int(record.label[3:])
        if segment[0] <= index <= segment[1]:
            residue_type = record.label[:3]
            color = RESIDUE_COLORS[residue_type]

            handle = ax.scatter(record[1], record[2], c=color)
            legend[residue_type] = handle

            if not nolabels:
                ax.annotate(f'{residue_type}{index + offset}', record[1:])
        
    # Show the legend when no labels are shown, or showlegend argument is passed
    if nolabels or showlegend:
        # Residues are added as they are encountered, not in sorted order
        sort_legend = {color: residue for color, residue in sorted(legend.items())}
        plt.legend(sort_legend.values(), sort_legend.keys())

    # Invert axes
    ax.invert_xaxis()
    ax.invert_yaxis()

    # Label axes
    atoms = record._fields[1:]
    ax.set_xlabel(AXIS_LABELS[atoms[0][0]])
    ax.set_ylabel(AXIS_LABELS[atoms[1][0]])

def plot3D(records: List[NamedTuple], *, nolabels: bool = False, showlegend: bool = False, offset: int, segment: Tuple[int], project: int, slices: int) -> None:
    
    # Get mins and maxes so all plots have the same xy coords
    maxes = [max(dim) for dim in list(zip(*records))[1:]]
    mins = [min(dim) for dim in list(zip(*records))[1:]]

    # Find intervals for bins
    projection_max = maxes.pop(project-1)
    projection_min = mins.pop(project-1)
    bin_width = (projection_max - projection_min) / slices
    cutoffs = [projection_min + i * bin_width for i in range(slices)]
    # Add 1 to not have to deal with < vs <= for last bin
    intervals = list(zip(cutoffs, cutoffs[1:] + [projection_max + 1]))

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
        for point in sorted(records, key=itemgetter(project)):
            residue_type = point[0][:3]
            index = int(point[0][3:])
            if low <= point[project] < high:
                # Remove projected index from list of indices used for plotting
                plot_point = list(point[:])
                plot_point.pop(project)
                sliced.append([f'{residue_type}{index+offset}'] + plot_point[1:])
        sliced_data.append(sliced)
    
    atoms = list(point._fields[1:])
    atoms.pop(project-1)
    ax.set_xlabel(AXIS_LABELS[atoms[0][0]])
    ax.set_ylabel(AXIS_LABELS[atoms[1][0]])

    tracker = Slices3D(ax, sliced_data, intervals)
    fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
    return tracker


class Slices3D:
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
                color = RESIDUE_COLORS[name[:3]]
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
                color = RESIDUE_COLORS[name[:3]]
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
import logging
from nightshift import constants
from typing import Dict, List, Tuple

import matplotlib.cm
from matplotlib.axes import Axes
from matplotlib.colors import is_color_like

AXIS_LABELS: Dict[str,str] = {
    'H': r'$^{1}$H (ppm)',
    'N': r'$^{15}$N (ppm)',
    'C': r'$^{13}$C (ppm)',
    }
Correlations = List[Tuple[int,str,Tuple[float]]]

def plot2D(
          ax: Axes, 
          correlations: Correlations,
          atoms: Tuple[str], 
          *,
          color: str = 'tab20',
          nolabels: bool = False,
          showlegend: bool = False,
          offset: int = 0) \
          -> Tuple[List,List]:
    legend = {}
    handles = []
    text = []
    residue_colors = get_residue_colors(color)
    
    for sequence_number, residue_type, chemical_shifts in correlations:
        color = residue_colors[residue_type]
        handle = ax.plot(*chemical_shifts, 'o', c=color)[0]
        handles.append(handle)
        legend[residue_type] = handle
        if not nolabels:
            text.append(ax.annotate(f'{residue_type}{sequence_number + offset}', chemical_shifts))
        
    # Show the legend when no labels are shown, or showlegend argument is passed
    if nolabels or showlegend:
        # Residues are added as they are encountered, not in sorted order
        sort_legend = {color: residue for residue, color in sorted(legend.items())}
        ax.legend(sort_legend.keys(), sort_legend.values())

    if not ax.xaxis_inverted():
        ax.invert_xaxis()
    if not ax.yaxis_inverted():
        ax.invert_yaxis()

    # Label axes
    ax.set_xlabel(AXIS_LABELS[atoms[0][0]])
    ax.set_ylabel(AXIS_LABELS[atoms[1][0]])
    
    return handles, text

def get_residue_colors(color: str) -> Dict:
    try:
        cmap = matplotlib.cm.get_cmap(color)
        color_vals = [cmap(i/20) for i in range(20)]
    except ValueError as err:
        if is_color_like(color):
            color_vals = [color]*20
        else:
            logging.warn(
                        f"color: '{color}' could not be interpreted as a color, setting" 
                        f"to default. If you wanted a colormap {str(err).partition(';')[-1]}"
                        )
            cmap = matplotlib.cm.get_cmap('tab20')
            color_vals = [cmap(i/20) for i in range(20)]
    return dict(zip(constants.ONE_LETTER_TO_THREE_LETTER.values(), color_vals))

class Slices3D:
    # Modified from: https://matplotlib.org/stable/gallery/event_handling/image_slices_viewer.html
    def __init__(self, ax, data, intervals, **kwargs):
        self.ax = ax
        self.data = data
        self.intervals = intervals
        self.slices = len(self.data)
        self.ind = self.slices//2 # start in the center
        self.text = []
        self.handles = []
        self.kwargs = kwargs
        self.update()

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
        self.clear_annotations()
        self.handles, self.text = plot2D(self.ax, self.data[self.ind], **self.kwargs)
        # Add ppm of center projected dimension
        self.text.append(self.ax.annotate(
                        f'{(sum(self.intervals[self.ind]))/2:.1f} ppm',
                        xy=(0.02,0.94), xycoords='axes fraction',
                        horizontalalignment='left'
                        ))
        self.ax.figure.canvas.draw()

    def clear_annotations(self):
        for text in self.text:
            text.remove()
            self.text = []
        for handle in self.handles:
            handle.set_data(None,None)
            self.handles = []

def plot3D(
          ax: Axes,
          correlations: Correlations,
          atoms: Tuple[str],
          *,
          color: str = 'tab20',
          nolabels: bool = False,
          showlegend: bool = False,
          offset: int = 0,
          project: int = 2,
          slices: int = 16) \
          -> Slices3D:
    
    # Get mins and maxes so all plots have the same xy coords
    # correlations formated [sequence_number, residue_type, (chemical_shifts)]
    shifts = list(zip(*(c[2] for c in correlations)))
    maxes = [max(dim) for dim in shifts]
    mins = [min(dim) for dim in shifts]

    # Find intervals for bins
    projection_max = maxes.pop(project)
    projection_min = mins.pop(project)
    bin_width = (projection_max - projection_min) / slices
    cutoffs = [projection_min + i * bin_width for i in range(slices)]
    # Add 1 to not have to deal with < vs <= for last bin
    intervals = list(zip(cutoffs, cutoffs[1:] + [projection_max + 1]))

    # Sets axes equal for all plots
    x_padding = 0.1 * (maxes[0] - mins[0])
    y_padding = 0.1 * (maxes[1] - mins[1])
    ax.set_xlim((mins[0]-x_padding, maxes[0]+x_padding))
    ax.set_ylim((mins[1]-y_padding, maxes[1]+y_padding))

    sliced_data = []
    for low, high in intervals:
        sliced = []
        # chemical shifts are index 2
        for sequence_number, residue_type, chemical_shifts in sorted(correlations, key=lambda x: x[2][project]):
            if low <= chemical_shifts[project] < high:
                # Remove projected index from list of indices used for plotting
                plot_point = list(chemical_shifts[:])
                plot_point.pop(project)
                sliced.append([sequence_number, residue_type, plot_point])
        sliced_data.append(sliced)
    
    atoms = atoms[:project] + atoms[project+1:]

    three_d = Slices3D(ax, sliced_data, intervals, atoms=atoms, color=color, nolabels=nolabels, showlegend=showlegend, offset=offset)
    ax.get_figure().canvas.mpl_connect('scroll_event', three_d.on_scroll)
    return three_d
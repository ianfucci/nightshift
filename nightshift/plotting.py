from typing import Dict, List, Tuple

from matplotlib.axes import Axes

# Palette of colors for each residue
RESIDUE_COLORS: Dict[str,str] = {
    'ALA':'maroon', 'ARG':'red', 'ASN':'pink', 'ASP':'brown',
    'CYS':'orange', 'GLU':'coral', 'GLN':'olive', 'GLY':'magenta',
    'HIS':'khaki', 'ILE':'purple', 'LEU':'green', 'MET':'navy',
    'LYS':'blue', 'PHE':'lime', 'PRO':'lightgreen', 'SER':'aquamarine',
    'THR':'cyan', 'TRP':'black', 'TYR':'grey', 'VAL':'yellow',
    }

AXIS_LABELS: Dict[str,str] = {'H': r'$^1$H (ppm)', 'N': r'$^{15}$N (ppm)', 'C': r'$^{13}$C (ppm)'}

def plot2D(ax: Axes, correlations: List, atoms: Tuple[str], 
          *, nolabels: bool, showlegend: bool, offset: int = 0) -> Tuple[List,List]:
    legend = {}
    handles = []
    text = []
    for sequence_number, residue_type, chemical_shifts in correlations:
        color = RESIDUE_COLORS[residue_type]
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

def plot3D(ax: Axes, correlations: List, atoms: Tuple[str], *, nolabels: bool = False,
          showlegend: bool = False, offset: int = 0, project: int, slices: int) -> None:
    
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
    
    atoms.pop(project)
    # ax.set_xlabel(AXIS_LABELS[atoms[0][0]])
    # ax.set_ylabel(AXIS_LABELS[atoms[1][0]])

    tracker = Slices3D(ax, sliced_data, intervals, atoms=atoms, nolabels=nolabels, showlegend=showlegend, offset=offset)
    ax.get_figure().canvas.mpl_connect('scroll_event', tracker.on_scroll)
    return tracker


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
        # self.handles, self.text = plot2D(self.ax, self.data[self.ind], **self.kwargs)
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
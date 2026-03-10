# -*- coding: utf-8 -*-
"""
miscellaneous utilities that do not belong anywhere else
"""
from __future__ import annotations
from mooonpy.tools.tables import ListListTable
from mooonpy.tools.file_utils import Path
import numpy as np
from typing import TYPE_CHECKING, Union, List, Optional
if TYPE_CHECKING:
    from matplotlib.axes import Axes

class ColorMap(ListListTable):
    """
    Read csv files of RGB values with a similar format to matplotlib.colormaps.
    root/tools/lookups contains some example custom colormaps, more can be added by adding files to this directory.

    **Currently Supported Colormaps**
        - coolwarm: blue to light to red diverging map (same as matplotlib?)
        - cooloverwarm: blue to grey to red diverging map
        - coswan: pink to grey to green diverging map
        - pyrolyze: green to grey to dark purple diverging map
        - redblack: red to dark grey linearish map
    """

    def __init__(self, name='coolwarm', low=0.0, high=1.0):
        """
        Read csv files of RGB vales into a ListListTable object.
        """
        super().__init__()
        thisfile = Path(__file__)
        parent = ListListTable.read_csv(thisfile.dir() / f'lookups/{name}.csv', True, True)
        self.__dict__.update(parent.__dict__)
        self.name = name
        self.low = low
        self.high = high

    def __call__(self, scalar, interp=False):
        """
        Return nearest color on 0-1 scale
        """

        if isinstance(scalar, float):
            scalar = (scalar - self.low) / (self.high - self.low)
            if scalar <= 0:
                row_ind = 0
            elif scalar >= 1:
                row_ind = -1
            else:
                row_ind = int(scalar * (len(self.grid) - 1))
            return tuple(self.grid[row_ind])
        elif isinstance(scalar, int):
            return tuple(self.grid[scalar])
        else:  # iterable
            out = []
            for scalar_ in scalar:
                scalar_ = float(scalar_ - self.low) / (self.high - self.low)
                if scalar_ <= 0:
                    row_ind = 0
                elif scalar_ >= 1:
                    row_ind = -1
                else:
                    row_ind = int(scalar_ * (len(self.grid) - 1))
                out.append(tuple(self.grid[row_ind]))
            return out


def inc_arrange(low, high, step):
    return np.arange(low, high + step, step)


def auto_colorbar(fig, axs, colormap, step, low=None, high=None, shrink=None, pad=None,
                  label='', orientation='vertical', kwargs=None):
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import ListedColormap, BoundaryNorm
    if isinstance(colormap, str):
        try:
            colormap = ColorMap(colormap)  # custom
        except:
            from matplotlib import colormaps as cm
            colormap = cm[colormap]  # mpl built in
            cmap = colormap

    if isinstance(colormap, ColorMap):
        if low is None:
            low = colormap.low
        if high is None:
            high = colormap.high
        ticks = inc_arrange(low, high, step)
        colors = colormap(ticks)
        cmap = ListedColormap(colors)
    else:
        if low is None:
            low = 0.0
        if high is None:
            high = 1.0
        ticks = inc_arrange(low, high, step)
        cmap = colormap

    bounds = np.arange(low - step / 2, high + step, step)
    norm = BoundaryNorm(bounds, cmap.N)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    if kwargs is None:
        kwargs = {}
    if shrink is not None:
        kwargs['shrink'] = shrink
    if pad is not None:
        kwargs['pad'] = pad
    cbar = fig.colorbar(sm, ax=axs, boundaries=bounds, ticks=ticks,
                        orientation=orientation, label=label, **kwargs)


def reset_axis_limits(
        ax: Axes,
        x: bool = True,
        y: bool = True,
        view: Optional[str] = None,
        style: Optional[Union[str, List[str]]] = None,
        xtick_spacing: Optional[float] = None,
        ytick_spacing: Optional[float] = None
) -> None:
    """
    Reset axis limits based on plotted line data with optional view presets and styling.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes object to modify
    x : bool
        Tighten x-axis limits to limit of data(default: True)
    y : bool
        Tighten y-axis limits to limit of data(default: True)
    view : str, optional
        Preset view configurations:
        - 'q1': First quadrant (x_min=0, y_min=0)
        - 'q2': Second quadrant (x_max=0, y_min=0)
        - 'q3': Third quadrant (x_max=0, y_max=0)
        - 'q4': Fourth quadrant (x_min=0, y_max=0)
        - '+x': Positive x (x_min=0, uses y parameter)
        - '-x': Negative x (x_max=0, uses y parameter)
        - '+y': Positive y (y_min=0, uses x parameter)
        - '-y': Negative y (y_max=0, uses x parameter)
    style : str or list of str, optional
        Styling options:
        - 'equal': Equal aspect ratio (ax.axis('equal'))
        - 'square': Square plot area (ax.axis('square'))
        - 'scaled': Scaled aspect ratio (ax.axis('scaled'))
        - 'tight': Tight layout (ax.axis('tight'))
        - 'auto': Auto scaling (ax.axis('auto'))
        - 'off': Turn off axis (ax.axis('off'))
        - 'on': Turn on axis (ax.axis('on'))
        - 'equalticks': Equal tick spacing with locked aspect ratio
        - 'box': Add box around plot (ax.set_frame_on(True))
        - 'nobox': Remove box (ax.set_frame_on(False))
        - 'grid': Add grid (ax.grid(True))
        - 'nogrid': Remove grid (ax.grid(False))
        Can pass a single string or list of strings for multiple styles
    xtick_spacing : float, optional
        Major tick spacing for x-axis (uses MultipleLocator)
    ytick_spacing : float, optional
        Major tick spacing for y-axis (uses MultipleLocator)
    """
    lines = ax.get_lines()
    collections = ax.collections  # For scatter plots
    patches = ax.patches

    if not lines and not collections and not patches:
        return  # No data to fit to

    # Get min/max x and y values from all lines and scatter plots
    x_data = []
    y_data = []

    # Get data from lines
    for line in lines:
        x_data.extend(line.get_xdata())
        y_data.extend(line.get_ydata())

    # Get data from scatter plots (PathCollections)
    for collection in collections:
        offsets = collection.get_offsets()
        if len(offsets) > 0:
            x_data.extend(offsets[:, 0])
            y_data.extend(offsets[:, 1])

    for patch in patches:
        # Get the bounding box of each patch
        bbox = patch.get_bbox()
        x_data.extend([bbox.x0, bbox.x1])
        y_data.extend([bbox.y0, bbox.y1])

    if not x_data or not y_data:
        return  # No valid data

    x_min = min(x_data)
    x_max = max(x_data)
    y_min = min(y_data)
    y_max = max(y_data)

    # Store current limits for view presets
    curr_xlim = ax.get_xlim()
    curr_ylim = ax.get_ylim()

    # Apply view presets
    if view:
        view = view.lower()

        # Set bounds to 0 based on view, only for the axes that are True
        if view in ['q1', 'q4', '+x']:
            x_min = 0
            if not x:
                x_max = curr_xlim[1]
        elif view in ['q2', 'q3', '-x']:
            x_max = 0
            if not x:
                x_min = curr_xlim[0]

        if view in ['q1', 'q2', '+y']:
            y_min = 0
            if not y:
                y_max = curr_ylim[1]
        elif view in ['q3', 'q4', '-y']:
            y_max = 0
            if not y:
                y_min = curr_ylim[0]

        # Set the axis limits to the data limits
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

    else:
        if x:
            ax.set_xlim(x_min, x_max)
        if y:
            ax.set_ylim(y_min, y_max)

    # Set tick spacing using MultipleLocator
    if xtick_spacing is not None or ytick_spacing is not None:
        from matplotlib.ticker import MultipleLocator
    if xtick_spacing is not None:
        ax.xaxis.set_major_locator(MultipleLocator(xtick_spacing))
    if ytick_spacing is not None:
        ax.yaxis.set_major_locator(MultipleLocator(ytick_spacing))

    # Apply styling
    if style:
        # Handle single string or list of strings
        if isinstance(style, str):
            style = [style]

        for s in style:
            s = s.lower()
            if s in ['equal', 'square', 'scaled', 'tight', 'auto', 'off', 'on']:
                ax.axis(s)
            elif s == 'equalticks':
                # Force matplotlib to update ticks
                ax.figure.canvas.draw()

                # Get the tick locations
                x_ticks = ax.get_xticks()
                y_ticks = ax.get_yticks()

                # Calculate tick spacing (assuming uniform spacing)
                if len(x_ticks) > 1 and len(y_ticks) > 1:
                    x_tick_spacing = x_ticks[1] - x_ticks[0]
                    y_tick_spacing = y_ticks[1] - y_ticks[0]

                    # Set aspect ratio so that tick spacing appears equal
                    # physical_spacing_x / physical_spacing_y = 1
                    # (x_tick_spacing / aspect) / y_tick_spacing = 1
                    # aspect = x_tick_spacing / y_tick_spacing
                    aspect_ratio = x_tick_spacing / y_tick_spacing
                    ax.set_aspect(aspect_ratio, adjustable='box')

            elif s == 'box':
                ax.set_frame_on(True)
            elif s == 'nobox':
                ax.set_frame_on(False)
            elif s == 'grid':
                ax.grid(True)
            elif s == 'nogrid':
                ax.grid(False)


def calculate_axis_ticks(data_min, data_max, ax=None, xaxis=True,
                         tick_range=(5, 9),
                         allowed_spacings=None, raise_error=True):
    """
    Calculate optimal axis limits and tick positions.

    Parameters
    ----------
    data_min : float
        Lower limit (usually 0)
    data_max : float
        Maximum data value to accommodate
    tick_range : tuple of int, optional
        (min_ticks, max_ticks), default (5, 9)
    allowed_spacings : list of float, optional
        Acceptable tick spacings, default [5, 10, 25, 50, 100, 200, 250]

    Returns
    -------
    ticks : ndarray
        Array of tick positions
    limit : float
        Upper axis limit
    spacing : float
        Spacing between ticks

    Raises
    ------
    ValueError
        If no valid tick configuration exists
    """
    if allowed_spacings is None:
        allowed_spacings = [0.05,0.1,0.25,0.5,1,5, 10, 15, 20, 25, 50, 100, 200, 250]

    min_ticks, max_ticks = tick_range

    # Find all valid combinations
    valid_configs = []

    for spacing in allowed_spacings:
        # Calculate how many ticks we'd need
        n_ticks = int(np.ceil((data_max - data_min) / spacing)) + 1

        # Check if within acceptable range
        if min_ticks <= n_ticks <= max_ticks:
            upper_limit = data_min + (n_ticks - 1) * spacing
            gap = upper_limit - data_max
            valid_configs.append({
                'spacing': spacing,
                'n_ticks': n_ticks,
                'limit': upper_limit,
                'gap': gap
            })

    if not valid_configs:
        # Calculate the extreme cases
        smallest_spacing = min(allowed_spacings)
        largest_spacing = max(allowed_spacings)

        # Smallest possible data_max: smallest spacing with min ticks
        min_data_max = data_min + (min_ticks - 1) * smallest_spacing

        # Largest possible data_max: largest spacing with max ticks
        max_data_max = data_min + (max_ticks - 1) * largest_spacing

        if raise_error:
            # Determine which constraint is violated
            if data_max > max_data_max:
                # Data range too LARGE
                raise ValueError(
                    f"No valid tick configuration found for data_min={data_min}, data_max={data_max} "
                    f"with tick_range={tick_range}.\n"
                    f"Data range too LARGE: maximum feasible data_max is {max_data_max:.4g}\n"
                    f"(using largest spacing={largest_spacing}, {max_ticks} ticks)"
                )
            elif data_max < min_data_max:
                # Data range too SMALL
                raise ValueError(
                    f"No valid tick configuration found for data_min={data_min}, data_max={data_max} "
                    f"with tick_range={tick_range}.\n"
                    f"Data range too SMALL: minimum feasible data_max is {min_data_max:.4g}\n"
                    f"(using smallest spacing={smallest_spacing}, {min_ticks} ticks)"
                )
            else:
                raise ValueError(
                    f"No valid tick configuration found for data_min={data_min}, data_max={data_max} "
                    f"with tick_range={tick_range}.\n"
                    f"Spacing infeasible:{allowed_spacings}\n"

                )
        else:
            return None, None, None
    # Sort by gap (ascending), then by n_ticks (ascending) for tiebreaking
    valid_configs.sort(key=lambda x: (x['gap'], x['n_ticks']))

    # Select best configuration
    best = valid_configs[0]

    # Generate tick array
    ticks = np.arange(data_min, best['limit'] + best['spacing'] / 2, best['spacing'])

    if ax is not None:
        if xaxis:
            ax.set(xticks=ticks, xlim=[data_min,best['limit']])
        else:
            ax.set(yticks=ticks, ylim=[data_min,best['limit']])
    return ticks, best['limit'], best['spacing']

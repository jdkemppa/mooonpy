# -*- coding: utf-8 -*-
import numpy as np
import warnings
from typing import Optional, Union

from mooonpy.tools.tables import ColTable
from mooonpy.tools.file_utils import Path
from mooonpy.tools.string_utils import _col_convert
# from ._files_io.read_logfile import readlog_basic

class Thermospace(ColTable):
    """
    Class to hold data found in LAMMPS logs.
    Data is organized into columns

    """

    def __init__(self, **kwargs):
        super(Thermospace, self).__init__(**kwargs)
        self.sections = {}

    def sect(self, sect_string: Optional[Union[str,range,list,int]]=None,priority='off') -> ColTable:
        if isinstance(sect_string, str):  # TODO
            if sect_string in self.sections:
                sections = [sect_string]
            else:
                pass  # Fancy splitting
        elif isinstance(sect_string, range):  # TODO
            sections = sect_string
        elif isinstance(sect_string, int):
            sections = [int(sect_string)]
        elif sect_string is None:
            sections = list(self.sections.keys())
        else:  # or allow custom ranges? also reverse slicing. look at slice object
            pass  # error message  # TODO

        out_table = ColTable(title=self.title, cornerlabel=self.cornerlabel)
        out_table.x_column = self.x_column
        if hasattr(self, 'default'):
            out_table.default = self.default

        for key, col in self.grid.items():
            sect_ranges = [self.sections[section] for section in sections]
            slices = []
            for ii, sect_range in enumerate(sect_ranges):
                if priority == 'first' and ii > 0:
                    slices.append(col[sect_range.start + sect_range.step:sect_range.stop:sect_range.step])
                elif priority == 'last' and ii > 0:
                    slices.append(col[sect_range.start:sect_range.stop - sect_range.step:sect_range.step])
                else: #  priority == 'off':
                    slices.append(col[sect_range.start:sect_range.stop:sect_range.step])
            out_table[key] = np.concatenate(slices)

        return out_table

    ## add merging options, and remove repeats method

    def __len__(self) -> Optional[int]:
        return self.shape()[0]

    @classmethod
    def basic_read(cls, file: Union[Path, str], silence_error_line: bool = False) -> 'Thermospace':
        return readlog_basic(file, silence_error_line=silence_error_line)

    def join_restart(self, restart, step_col='Step',restart_step_ind=0, this_step_ind=None):
        """
        Appends data from a restarted thermospace to the current one.
        Uses (step_col, restart_step_ind) as the start of the appended section,
        and finds the last matching (step_col, this_step_ind) and overwrites after that

        step_col is the column used for indexing, defaults to 'Step' but 'v_Time' may also be used

        restart_step_ind = 0 uses 0th index step, and so on
        this_step_ind = None uses last matching value
        """
        this_x = self[step_col]
        restart_x = restart[step_col]
        if this_step_ind is None:
            matches = np.argwhere(np.equal(this_x, restart_x[restart_step_ind]))
            if len(matches) == 0:
                this_step_ind = len(this_x)  # no slice
            else:
                this_step_ind = matches[-1][0] # last match

        for key in self.grid.keys():
            if key in restart.grid:
                self.grid[key] = np.concatenate((self.grid[key][:this_step_ind], restart.grid[key][restart_step_ind:]))
            else:
                self.grid[key] = np.concatenate((self.grid[key][:this_step_ind],np.full(len(restart),np.nan)[restart_step_ind:]))

        last_sect = None
        for sect, range_ in self.sections.items():
            if this_step_ind in range_:
                last_sect = sect
                # Continues to last match
        if last_sect is None:
            self.sections.update(restart.sections)
        else:
            for sect, range_ in restart.sections.items():
                if restart_step_ind in range_:
                    self.sections[last_sect] = range(self.sections[last_sect].start, range_.stop+this_step_ind)
                else:
                    self.sections[last_sect+sect-1] = range(range_.start+this_step_ind, range_.stop+this_step_ind)
        ## Test this ^
        # !!!

## This should be refactored into _files_io but imports are being weird
def readlog_basic(file: [Path, str], silence_error_line: bool = False) -> Thermospace:
    """
    Read a single log file into a Thermospace object.
    Only reads thermo table, no timing or variable information.

    :param file: path to a log file
    :type file: [Path,str]
    :param silence_error_line: silences error line and warnings if True (default False)
    :type silence_error_line: bool
    :return: Thermospace object
    :rtype: Thermospace

    :Example:
        >>> import mooonpy
        >>> file = mooonpy.Path('somepath.log.lammps')
        >>> MyLog = mooonpy.readlog_basic(file)
        >>> MyLog.csv(file.new_ext('.csv'))

    .. seealso:: :class:`thermospace.Thermospace`, :class:`tables.ColTable`
    .. warning:: Files with no thermo data, or incorrect format, will raise a warning, then return an empty Thermospace object.
    .. note:: This function is capable reading logs with pauses for XRD, and changes to
        'LAMMMPS thermo_style'_. Columns with missing data from style changes are padded with
        np.nan values in a float array. Complete columns that are intergers are converted to int arrays.

    .. todo::
        - Unit tests for failure modes and column changes

    .. _LAMMMPS thermo_style: https://docs.lammps.org/thermo_style.html

    """
    ## variables to return in thermospace
    file = Path(file)
    columns = {}
    sections = {}
    ## Setup for internal variables
    has_nan = set()
    keywords = []  # dummy
    rowindex = 0  # starts at index 0
    sectionID = 0
    startrow = 0  # dummy
    data_flag = False
    interrupt_flag = False
    header_flag = True

    if not file:
        raise Exception(f'File {file} not found')
    with file.open('r') as f:
        for line in f:
            line = line.strip()
            ## Exit conditions and switch cases
            if not line:
                continue
            elif not line[0].isdigit():  # single check is cheaper than 5
                if line.startswith('WARNING'):
                    continue # add warning logger to advanced version
                elif line.startswith('ERROR'):
                    if not silence_error_line:
                        print('File {:} contains Error line, exiting read'.format(file))
                    break
                elif 'Sending Ctrl-C to processes as requested' in line:
                    if not silence_error_line:
                        print('File {:} contains Ctrl-C exit, exiting read'.format(file))
                    break
                elif line.startswith('Per MPI'):
                    data_flag = True
                    header_flag = True
                    continue
                elif line.startswith('Loop time of'):
                    data_flag = False
                    sections[sectionID] = range(startrow, rowindex)  ## check these indexes
                    for missing in set(columns.keys()).difference(set(keywords)):
                        columns[missing] += [None] * (rowindex - startrow)
                        has_nan.add(missing)
                    continue
                elif line == '-----':  ## for XRD sims, not sure what else
                    interrupt_flag = not interrupt_flag  ## Toggle reading
                    continue

            ## Read thermo section
            if data_flag and not interrupt_flag:
                splits = line.split()
                if header_flag:
                    header_flag = False
                    keywords = splits
                    sectionID += 1
                    startrow = rowindex
                    for key in keywords:
                        if key not in columns:
                            columns[key] = [None] * (rowindex)  ## init values for new columns, [] if rowcount is -1
                            if rowindex != 0:
                                has_nan.add(key)
                else:  # table body
                    if len(keywords) != len(splits):
                        if not silence_error_line:
                            print('File {:} ends unexpectedly skipping last line'.format(file))
                        break
                    for k, v in zip(keywords, splits):
                        columns[k].append(v)  ## no string to float conversion, handled by numpy conversion later
                    rowindex += 1  ## after in case anything fails
            ## End thermo block
        ## End read loop
    ## End with statement
    sections[sectionID] = range(startrow, rowindex + 1)
    # ^ add section for last successful row, and increase final index by 1

    for missing in set(columns.keys()).difference(set(keywords)):  # add nan to missing columns
        columns[missing] += [None] * (rowindex - startrow)
        has_nan.add(missing)

    for key, col in columns.items():  # convert string lists to array
        nan = bool(key in has_nan)
        # print(key, nan)
        columns[key] = _col_convert(col, nan)
    if len(columns) == 0 and not silence_error_line:
        warnings.warn(f'File {file} Contains no thermo data.')
    out = Thermospace()  ## may change with init?
    out.grid = columns
    out.title = file
    out.sections = sections
    return out

if __name__ == '__main__':
    this_x = np.array([0,1,2,3])
    restart_x = np.array([3,4,5,6,7])
    restart_step_ind =0
    matches = np.argwhere(np.equal(this_x,restart_x[restart_step_ind]))
    if len(matches) == 0:
        this_step_ind = len(this_x) # no slice
    else:
        this_step_ind = matches[-1][0]

    new = np.concatenate((this_x[:this_step_ind], restart_x[restart_step_ind:]))
    print(new)
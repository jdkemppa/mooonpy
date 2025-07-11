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
    Data is organized into

    """

    def __init__(self, **kwargs):
        super(Thermospace, self).__init__(**kwargs)
        self.sections = {}

    def sect(self, sect_string: Union[str,range,list,int]) -> ColTable:
        if isinstance(sect_string, str):  # TODO
            if sect_string in self.sections:
                sections = [sect_string]
            else:
                pass  # Fancy splitting
        elif isinstance(sect_string, range):  # TODO
            pass
        elif isinstance(sect_string, int):
            sections = [int(sect_string)]
        else:  # or allow custom ranges? also reverse slicing. look at slice object
            pass  # error message  # TODO

        out_table = ColTable(title=self.title, cornerlabel=self.cornerlabel)
        out_table.x_column = self.x_column
        if hasattr(self, 'default'):
            out_table.default = self.default

        for key, col in self.grid.items():
            sect_ranges = [self.sections[section] for section in sections]
            column = np.concatenate(
                [self[key][sect_range.start:sect_range.stop:sect_range.step] for sect_range in sect_ranges])
            out_table[key] = column

        return out_table

    ## add merging options, and remove repeats method

    def __len__(self) -> Optional[int]:
        return self.shape()[0]

    @classmethod
    def basic_read(cls, file: Union[Path, str], silence_error_line: bool = False) -> 'Thermospace':
        return readlog_basic(file, silence_error_line=silence_error_line)

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
                # if line.startswith('WARNING'):
                #     continue
                if line.startswith('ERROR'):
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
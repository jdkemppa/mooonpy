# -*- coding: utf-8 -*-
import numpy as np
import warnings
from typing import Optional, Union

from mooonpy.tools.tables import ColTable
from mooonpy.tools.file_utils import Path
from mooonpy.tools.string_utils import _col_convert


class Thermospace(ColTable):
    """
    Class to hold data found in LAMMPS logs.
    Data is organized into columns
    """

    def __init__(self, filename=None, **kwargs):
        super(Thermospace, self).__init__(**kwargs)
        self.sections = {}

        if filename is not None:  # populate from file
            readlog_basic(filename, out=self)

    def sect(self, sect_string: Optional[Union[str, range, list, int, np.ndarray]] = None, priority='off') -> ColTable:
        """
        Extract a subset of data by section ID(s).

        :param sect_string: Section selector.
            int for single section, list/array of ints for multiple,
            range for a range of IDs, None for all sections.
            str supports comma-separated values and colon ranges with stride
            (inclusive, 1-indexed, negative indexing):
            ``"1:3"`` ``"6:"`` ``"-3:"`` ``"::-1"`` ``"1:5:2"`` ``"1,3,-1"``
        :param priority: 'first' drops the first row of non-first sections,
            'last' drops the last row of non-first sections, 'off' keeps all rows.
        """
        priority = priority.lower()
        if isinstance(sect_string, (int, np.integer)):
            sections = [int(sect_string)]
        elif isinstance(sect_string, np.ndarray):
            sections = sect_string.ravel().astype(int).tolist()
        elif isinstance(sect_string, (list, range)):
            sections = list(sect_string)
        elif isinstance(sect_string, str):
            sections = self._parse_sect_string(sect_string)
        elif sect_string is None:
            sections = list(self.sections.keys())
        else:
            raise TypeError(f'sect() expected int, list, range, str, array, or None; got {type(sect_string).__name__}')

        out_table = ColTable(title=self.title, cornerlabel=self.cornerlabel)
        out_table.x_column = self.x_column
        if hasattr(self, 'default'):
            out_table.default = self.default

        sect_ranges = []
        for section in sections:
            sect_range = self.sections.get(section)
            if sect_range is not None:
                sect_ranges.append(sect_range)

        if sect_ranges:
            for key, col in self.grid.items():
                slices = []
                for ii, sect_range in enumerate(sect_ranges):
                    start = sect_range.start
                    stop = sect_range.stop
                    if ii > 0 and priority == 'first':
                        start += 1  # drop first row of continuation sections
                    elif ii > 0 and priority == 'last':
                        stop -= 1  # drop last row of continuation sections
                    slices.append(col[start:stop])
                out_table[key] = np.concatenate(slices)

        return out_table

    def _parse_sect_string(self, sect_string: str) -> list:
        """
        Parse section string into list of section IDs.
        Comma-separated terms, each term is ``start:stop:step`` (any part optional).
        Endpoints are inclusive, 1-indexed, negative indices count from end.

        Examples: ``'6:'``, ``'-1'``, ``'-3:'``, ``'::-1'``, ``'1:5:2'``, ``'1,3,-1'``
        """
        max_sect = max(self.sections.keys()) if self.sections else 0

        def resolve(val, default):
            """Resolve a single index: parse int, apply negative indexing, or use default."""
            val = val.strip()
            if not val:
                return default
            n = int(val)
            if n < 0:
                return max_sect + 1 + n  # -1 → max_sect, -2 → max_sect-1
            return n

        sections = []
        for part in sect_string.split(','):
            part = part.strip()
            if ':' in part:
                pieces = part.split(':')
                step = resolve(pieces[2], 1) if len(pieces) > 2 else 1
                if step > 0:
                    start = resolve(pieces[0], 1)
                    stop = resolve(pieces[1], max_sect)
                else:
                    start = resolve(pieces[0], max_sect)
                    stop = resolve(pieces[1], 1)
                # inclusive endpoints
                sections.extend(range(start, stop + (1 if step > 0 else -1), step))
            else:
                sections.append(resolve(part, None))
        return sections

    ## add merging options, and remove repeats method

    def __len__(self) -> Optional[int]:
        return self.shape()[0]

    # @classmethod
    # def basic_read(cls, file: Union[Path, str], silence_error_line: bool = False) -> 'Thermospace':
    #     return readlog_basic(file, silence_error_line=silence_error_line)

    def join_restart(self, restart, step_col='Step', restart_step_ind=0, this_step_ind=None):
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
                this_step_ind = matches[-1][0]  # last match

        for key in self.grid.keys():
            if key in restart.grid:
                self.grid[key] = np.concatenate((self.grid[key][:this_step_ind], restart.grid[key][restart_step_ind:]))
            else:
                self.grid[key] = np.concatenate(
                    (self.grid[key][:this_step_ind], np.full(len(restart), np.nan)[restart_step_ind:]))

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
                    self.sections[last_sect] = range(self.sections[last_sect].start, range_.stop + this_step_ind)
                else:
                    self.sections[last_sect + sect - 1] = range(range_.start + this_step_ind,
                                                                range_.stop + this_step_ind)
        ## Test this ^
        # !!!


## This should be refactored into _files_io but imports are being weird
# from ._files_io.read_logfile import readlog_basic # put back later after more reads
def readlog_basic(file: [Path, str], silence_error_line: bool = False, out: Thermospace = None) -> Thermospace:
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
            elif not (line[0].isdigit() or (line[0] == '-' and line[1].isdigit())):  # data rows start with digit or negative sign
                if line.startswith('WARNING'):
                    continue  # add warning logger to advanced version
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
    if sectionID not in sections:  # close last section if not already closed by "Loop time of"
        sections[sectionID] = range(startrow, rowindex)
        for missing in set(columns.keys()).difference(set(keywords)):  # add nan to missing columns
            columns[missing] += [None] * (rowindex - startrow)
            has_nan.add(missing)

    for key, col in columns.items():  # convert string lists to array
        nan = bool(key in has_nan)
        # print(key, nan)
        columns[key] = _col_convert(col, nan)
    if len(columns) == 0 and not silence_error_line:
        warnings.warn(f'File {file} Contains no thermo data.')
    if not isinstance(out, Thermospace) or out is None:
        out = Thermospace()  ## may change with init?
    out.grid = columns
    out.title = file
    out.sections = sections
    return out


if __name__ == '__main__':
    this_x = np.array([0, 1, 2, 3])
    restart_x = np.array([3, 4, 5, 6, 7])
    restart_step_ind = 0
    matches = np.argwhere(np.equal(this_x, restart_x[restart_step_ind]))
    if len(matches) == 0:
        this_step_ind = len(this_x)  # no slice
    else:
        this_step_ind = matches[-1][0]

    new = np.concatenate((this_x[:this_step_ind], restart_x[restart_step_ind:]))
    print(new)

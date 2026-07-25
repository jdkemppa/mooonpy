# -*- coding: utf-8 -*-
import numpy as np
from typing import Optional, Union

from mooonpy.thermospace import _files_io
from mooonpy.tools.tables import ColTable
from mooonpy import Path


class Thermospace(ColTable):
    """
    Class to hold data found in LAMMPS logs.
    Data is organized into columns
    """

    def __init__(self, filename=None, **kwargs):
        silent = kwargs.pop('silence_error_line', False)
        super(Thermospace, self).__init__(**kwargs)
        self.sections = {}

        if filename is not None:  # populate from file
            readlog_basic(filename, out=self,silence_error_line=silent)

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


# from ._files_io.read_logfile import readlog_basic # put back later after more reads
def readlog_basic(file: Path | str, silence_error_line: bool = False, out: Thermospace = None) -> Thermospace:
    from mooonpy.thermospace import read_lammps_logfile
    print("'readlog_basic' has been renamed read_lammps_logfile and can now "
          "be called from mooonpy.thermospace import read_lammps_logfile.")
    return read_lammps_logfile(file, silence_error_line, out)


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

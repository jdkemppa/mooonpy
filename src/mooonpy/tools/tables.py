# -*- coding: utf-8 -*-
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from typing import Optional, Tuple, List, Dict

from .string_utils import _col_convert, string2digit
from .math_utils import aggregate_fun
from .file_utils import Path


class Table(object):
    """
    The Table object takes several form-factors of the underlying data structure
    implemented in its subclasses. The main class is a dummy for storage but
    handles IO to csv or text. Methods are developed so access can be aliased
    for any of the subclasses. Or the .grid attribute can be used for direct
    access on specific subclasses.

    def shape(self) -> Tuple[Optional[int], Optional[int]]:
        '''
        Size of the grid attribute, similar to numpy's .shape attribute.
        If one item in the tuple is None, the grid is asymmetric or does
        not match the labels and shape is meaningless.
        '''
        return tuple((0,0)) # makes code happy
    def headers(self) -> Optional[List[str]]:
        '''
        Return a list of keywords or headers of the table
        '''
        pass
    def rowlabels(self) -> Optional[List[str]]:
        '''
        Returns a list of row labels of the table
        '''
        pass
    def col(self,key):
        '''
        Returns column from index of 'key' in headers, or if not found, index of the header
        if that fails, a column of appropriate length is created with .default
        '''
        pass
    def row(self,row_key):
        '''
        Returns a row from the index of 'row_key' in rows, or if not found is used as the index
        This has O(N) cost if using 'rowlabel' as a list, or O(1) cost without the list
        '''
        pass
    def rowcol(self,row_key,key):
        pass

    """

    def __init__(self, title=None, cornerlabel=None, default=...):
        self.title = title
        self.cornerlabel = cornerlabel
        if default is not ...:  # most likely default is 0 or None so Ellipsis is used to skip the attribute creation
            self.default = default
        self.grid = None  ## this data structure changes for each subclass

        self.summary = {}
        # These are output formatting variables, currently not kwarg'd into existence
        self.delim = '\t'
        self.cell = '{:<8}'
        self.float_cell = '{:< 8.3f}'
        self.spacer = None
        self.header = True

        self.x_column = None  # used for .plot() method

    def col_fun(self, fun_name: str) -> Dict | List:
        """
        Compute aggregate functions on each column.

        The result is added to the summary dict as a 'fun_name':dict(key:value,...)

        pair if keywords are defined, or a list(value,...) if not.

        note that the summary dict may be converted to a RowTable

        RowTable does not exist yet

        The data is also returned.

        """
        keys = self.headers()
        if keys:
            out = dict()
            for key in keys:
                out[key] = aggregate_fun(fun_name, np.array(self.col(key)))
        else:
            out = list()
            for col_index in range(self.shape()[1]):
                out.append(aggregate_fun(fun_name, np.array(self.col(col_index))))
        self.summary[fun_name] = out
        return out

    def __str__(self):
        type_ = self.__class__.__name__
        shape = self.shape()
        title = self.title
        if self.title:
            return f'{type_} of size {shape} and title "{title}"'
        else:
            return f'{type_} of size {shape}'

    def console(self, **kwargs):
        """
        Convert table to string representation
        This is ok to use in console for small tables, but will struggle for file IO with large tables
        """
        self.__dict__.update(**kwargs)  # update with formatting kwargs
        for line in self.yield_line():
            print(line)

    def csv(self, file: Path | str, append=False, skip_header=False):
        """
        Write table to csv file
        """
        file = Path(file)
        self.delim = ','  # override
        self.spacer = None
        if file and append:
            mode = 'a'
        else:
            mode = 'w'
        with open(file, mode) as out_file:
            for line in self.yield_line():
                if skip_header: skip_header = False; continue
                out_file.write(line + '\n')

    def plot(self, **kwargs):
        self.__dict__.update(kwargs)  # update with formatting kwargs
        fig, axs = plt.subplots()
        headers = self.headers()

        if headers:
            if self.x_column is not None:
                x_column = self.x_column
                x = self.col(x_column)
                for header in headers:
                    if header == x_column: continue
                    axs.plot(x, self.col(header), label=header)
                axs.set_xlabel(x_column)
                axs.legend()
            else:
                for header in headers:
                    axs.plot(self.col(header), label=header)
                axs.legend()
        else:
            for col_index in range(self.shape()[1]):
                axs.plot(self.col(col_index))

        if self.title:
            axs.title.set_text(self.title)
        return fig, axs

    def yield_line(self):
        """
        Create one line of a table with a specified format
        Will probably only work for rectangular tables
        """
        left_flag = bool(self.rowlabels())  # skip for None or empty
        top_flag = bool(self.headers()) and self.header
        shape = self.shape()
        cell = self.cell  # keep self if changes outside?
        float_cell = self.float_cell
        delim = self.delim
        if self.spacer is not None:
            if left_flag:
                width = shape[1] + 1
            else:
                width = shape[1]
            space_line = self.spacer * width
            yield space_line  # 1st line
        else:
            space_line = ''

        # if self.title is not None:
        #     yield self.title
        #     if space_line: yield space_line

        line = ''
        if left_flag and top_flag:
            if self.cornerlabel:
                line += cell.format(self.cornerlabel) + delim
            else:
                line += cell.format('') + delim  # empty cell

        if top_flag:
            for key in self.headers():
                line += cell.format(key) + delim
            yield line

        if space_line: yield space_line

        for row_index in range(shape[0]):
            line = ''
            if left_flag:
                line += cell.format(self.rowlabels()[row_index]) + delim
            for col_index in range(shape[1]):
                if isinstance(self, ListListTable): # weird indexing error with rowcol here so it was hardcoded out
                    value = self.grid[row_index][col_index]
                else:
                    value = self.rowcol(row_index, col_index)
                if value is None:
                    line += delim
                elif isinstance(value, float):
                    line += float_cell.format(value) + delim
                else:
                    line += cell.format(value) + delim
            yield line

        if space_line: yield space_line

    def copy(self):
        return deepcopy(self)


class ArrayTable(Table):
    def __init__(self, from_numpy=None, shape=None, collabels=None, rows=None, title=None, cornerlabel=None,
                 default=...):
        super(ArrayTable, self).__init__(title, cornerlabel,
                                         default)  # dict creation is bad in function def, but it's not used for anything so probably fine?
        if from_numpy is not None:
            self.grid = np.array(from_numpy)
        elif shape is not None:
            if hasattr(self, 'default'):
                self.grid = np.full(shape, self.default)  # assume tuple shape
            else:
                self.grid = np.empty(shape)  # assume tuple
        else:
            self.grid = np.empty(0)

        if collabels is not None:
            self.keywords = collabels
        else:
            self.keywords = None
        if rows is not None:
            self.rows = rows
        else:
            self.rows = None

    def shape(self) -> Tuple[Optional[int], Optional[int]]:
        return self.grid.shape  # use numpy def

    def headers(self) -> Optional[List[str]]:
        if self.keywords is not None:
            return list(self.keywords)
        else:
            return None

    def rowlabels(self) -> Optional[List[str]]:
        if self.rows is not None:
            return list(self.rows)
        else:
            return None

    def col(self, key):
        try:
            if self.keywords is not None and key in self.keywords:
                index = list(self.keywords).index(key)
                return self.grid[:, index]
            else:
                return self.grid[:, int(key)]
        except:
            if hasattr(self, 'default'):
                return np.full(self.shape()[0], self.default)  # return array of default values
            else:
                raise Exception(f'ERROR: ColTable does not have {key} as a column label, index or a default attribute')

    def row(self, row_key):
        try:
            if self.rows is not None and row_key in self.rows:
                index = list(self.rows).index(row_key)
                return self.grid[:, index]
            else:
                return self.grid[:, int(row_key)]
        except:
            if hasattr(self, 'default'):
                return np.full(self.shape()[1], self.default)  # return array of default values
            else:
                raise Exception(f'ERROR: ColTable does not have {row_key} as a row label, index or a default attribute')

    def rowcol(self, row_key, key):
        if self.rows is None:
            row_index = int(row_key)
        elif row_key in self.rows:
            row_index = self.rows.index(row_key)
        else:
            row_index = int(row_key)

        if self.keywords is None:
            col_index = int(key)
        elif key in self.keywords:
            col_index = self.keywords.index(key)
        else:
            col_index = int(key)

        try:
            value = self.grid[row_index, col_index]  # numpy indexing
            return value
        except:
            if hasattr(self, 'default'):
                return self.default
        raise Exception(f'ERROR: ArrayTable row {row_key} and col {key} could not find a match or default')


class ListListTable(Table):
    def __init__(self, from_listlist=None, shape=None, collabels=None, rows=None, title=None, cornerlabel=None,
                 default=...):
        super(ListListTable, self).__init__(title, cornerlabel, default)
        if from_listlist is not None:
            self.grid = from_listlist
        elif shape is not None:
            if hasattr(self, 'default'):
                # self.grid = [[self.default] * shape[1]] * shape[0]
                self.grid = []
                for rr in range(shape[0]):
                    self.grid.append([])
                    for cc in range(shape[1]):
                        self.grid[rr].append(default)
            else:
                # self.grid = [[None] * shape[1]] * shape[0]  # empty does not make sense in this case, using None as default
                self.grid = []
                for rr in range(shape[0]):
                    self.grid.append([])
                    for cc in range(shape[1]):
                        self.grid[rr].append(None)
        else:
            self.grid = []

        if collabels is not None:
            self.keywords = collabels
        else:
            self.keywords = None
        if rows is not None:
            self.rows = rows
        else:
            self.rows = None

    def append(self, row):
        """Append to ehd of grid"""
        self.grid.append(row)

    def shape(self) -> Tuple[Optional[int], Optional[int]]:
        n_rows = len(self.grid)
        set_col = set([len(row) for row in self.grid])
        if len(set_col) == 0:
            n_cols = 0
        elif len(set_col) == 1:
            n_cols = set_col.pop()  # grab common length
        else:
            return tuple((n_rows, None))  # uncommon lengths exit here
        if self.keywords is not None:
            if len(self.keywords) != n_cols:
                n_cols = None  # asymmetry
        return tuple((n_rows, n_cols))

    def headers(self) -> Optional[List[str]]:
        if self.keywords is not None:
            return list(self.keywords)
        else:
            return None

    def rowlabels(self) -> Optional[List[str]]:
        if self.rows is not None:
            return list(self.rows)
        else:
            return None

    def col(self, key):
        try:
            if self.keywords is not None and key in self.keywords:
                index = list(self.keywords).index(key)
                return [row[index] for row in self.grid]
            else:
                index = int(key)
                return [row[index] for row in self.grid]
        except:
            if hasattr(self, 'default'):
                return np.full(self.shape()[0], self.default)  # return array of default values
            else:
                raise Exception(f'ERROR: ListListTable does not have {key} as a column label or a default attribute')

    def row(self, row_key):
        try:
            if self.rows is not None and row_key in self.rows:
                index = list(self.rows).index(row_key)
                return self.grid[index]
            else:
                return self.grid[int(row_key)]
        except:
            if hasattr(self, 'default'):
                return [self.default] * self.shape()[1]  # return list of default values
            else:
                raise Exception(
                    f'ERROR: ListListTable does not have {row_key} as a row label, index or a default attribute')

    def rowcol(self, row_key, key):
        if self.rows is None:
            row_index = int(row_key)
        elif row_key in self.rows:
            row_index = self.rows.index(row_key)
        else:
            row_index = int(row_key)

        if self.keywords is None:
            col_index = int(key)
        elif key in self.keywords:
            col_index = self.keywords.index(key)
        else:
            col_index = int(key)

        try:
            value = self.grid[row_index][col_index]
            return value
        except:
            if hasattr(self, 'default'):
                return self.default
        raise Exception(f'ERROR: ListListTable row {row_key} and col {key} could not find a match or default')

    @classmethod
    def read_csv(cls, file: Path | str, header=True, rowlabels=False):
        """
        Read csv file into columns
        """
        file = Path(file)
        file_obj = file.open(mode='r')
        row_list = []
        columns = None
        corner = None
        rows = []

        try:
            for ii, line in enumerate(file_obj):
                line = line.strip()
                splits = line.split(',')

                if ii == 0:
                    if header:
                        if rowlabels:
                            corner = splits[0].strip()
                            columns = [key.strip() for key in splits[1:]]
                            rows = []
                            if columns[-1] == '':
                                columns.pop(-1)  # remove trailing comma
                        else:
                            rows = None
                        continue # no data

                this_row = []
                if rowlabels:
                    rows.append(splits[0].strip())

                    for value in splits[1:]:
                        this_row.append(string2digit(value.strip()))
                else:
                    for value in splits:
                        this_row.append(string2digit(value.strip()))
                row_list.append(this_row)
        except:
            raise Exception(f'ERROR: CSV line {ii} :{line} could not read from {file}')

        return ListListTable(from_listlist=row_list, cornerlabel=corner, rows=rows,collabels=columns)

class ColTable(Table):
    def __init__(self, from_dict=None, rows=None, title=None, cornerlabel=None, default=...):
        super(ColTable, self).__init__(title, cornerlabel, default)
        if from_dict is not None:
            self.grid = dict(from_dict)  # works for constructors as well? Test this
        else:
            self.grid = {}
        self.rows = rows

    def __setitem__(self, key, value):
        self.grid[key] = value

    def __getitem__(self, key):  # alias
        return self.col(key)

    def __delitem__(self, key):
        del self.grid[key]

    def __contains__(self, key):
        return key in self.grid

    def shape(self) -> Tuple[Optional[int], Optional[int]]:
        n_cols = len(self.grid)  # will always exist >= 0
        if self.rows is not None:
            n_rows = len(self.rows)  # compare to rowlabel
        else:
            n_rows = None  # temporary
        for key, col in self.grid.items():
            if n_rows is None:
                n_rows = len(col)  # 1st col is used for compare
            elif len(col) != n_rows:
                n_rows = None
                break
            # else is matching
        return tuple((n_rows, n_cols))

    def headers(self) -> Optional[List[str]]:
        return list(self.grid.keys())  # numpy?

    def rowlabels(self) -> Optional[List[str]]:
        if self.rows is not None:
            return list(self.rows)  # for parity with other form factors
        else:
            return None

    def row(self, row_key):
        try:
            if self.rows is None:
                index = int(row_key)  # if lookup not defined, use index in columns
            elif row_key in self.rows:
                row_list = list(self.rows)  # might be numpy
                index = row_list.index(row_key)
            else:
                index = int(row_key)
            out_row = []
            for key, col in self.grid.items():
                out_row.append(col[index])
            return out_row
        except:
            if hasattr(self, 'default'):
                return self.default
        raise Exception(f'ERROR: ColTable row {row_key} could not find a match or default')

    def col(self, key):
        try:
            if key in self.grid:
                return self.grid[key]
            elif hasattr(self, 'default'):
                return np.full(self.shape()[0], self.default)  # return array of default values
        except:
            pass
        raise Exception(f'ERROR: ColTable does not have {key} as a column label or a default attribute')

    def rowcol(self, row_key, key):

        try:
            if key not in self.grid:
                key = list(self.grid.keys())[key]  # use as index
            if self.rows is None:
                value = self.grid[key][int(row_key)]  # assume int index
            elif row_key in self.rows:
                index = self.rows.index(row_key)
                value = self.grid[key][index]
            else:
                value = self.grid[key][int(row_key)]  # try int, or it will crash
            return value
        except:
            if hasattr(self, 'default'):
                return self.default
        raise Exception(f'ERROR: ColTable row {row_key} and col {key} could not find a match or default')

    @classmethod
    def read_csv(cls, file: Path | str, header=True, rowlabels=False, keep_list=None):
        """
        Read csv file into columns
        """
        file = Path(file)
        file_obj = file.open(mode='r')
        columns = {}
        corner = None
        rows = None
        if keep_list is None: keep_list = []

        try:
            for ii, line in enumerate(file_obj):
                line = line.strip()
                splits = line.split(',')

                if ii == 0:
                    if header:
                        if rowlabels:
                            corner = splits[0].strip()
                            keys = [key.strip() for key in splits[1:]]
                            rows = []
                        else:
                            keys = [key.strip() for key in splits]
                    else:
                        if rowlabels:
                            keys = list(range(len(splits) - 1))
                            rows = []
                        else:
                            keys = list(range(len(splits)))
                    if keys[-1] == '': keys.pop(-1) # remove trailing comma
                    columns = {key: list() for key in keys}
                else:
                    if rowlabels:
                        rows.append(splits[0].strip())
                        for key, value in zip(keys, splits[1:]):
                            columns[key].append(value.strip())
                    else:
                        for key, value in zip(keys, splits):
                            columns[key].append(value.strip())
        except:
            raise Exception(f'ERROR: CSV line {ii} :{line} could not read from {file}')

        for ii, (key, vector) in enumerate(columns.items()):
            if ii in keep_list or key in keep_list: continue
            try:
                columns[key] = _col_convert(vector)  # string to float converter
            except:
                pass  # keep as list of strings

        return ColTable(from_dict=columns, cornerlabel=corner, rows=rows)


# %%
if __name__ == "__main__":
    coltab = ColTable(from_dict={'Step': [0, 1, 2], 'Temp': [300, 475, 600]})
    print(coltab.shape())
    print(coltab.headers())
    print(coltab['Step'])
    print(coltab.col('Step'))  # alias for __getitem__
    print(coltab.row(0))  # index access
    print(coltab.rowlabels())  # empty
    coltab.rows = ['start', 'mid', 'end']
    print(coltab.rowlabels())
    print(coltab.row('mid'))  # .index in rows lookup
    print(coltab.row(0))  # index access
    # print(coltab.row(3))  # bad index access
    print(coltab.rowcol('end', 'Temp'))
    print(coltab.rowcol(2, 1))

    print()
    emptarr = ArrayTable(shape=(3, 2))
    print(emptarr.headers())

    print()
    arrtab = ArrayTable([[0, 300], [1, 475], [2, 600]])
    arrtab.col_fun('sum')
    arrtab.col_fun('sum0')
    arrtab.col_fun('avg')
    arrtab.col_fun('avg0')
    print(arrtab.shape())
    arrtab.keywords = ['Step', 'Temp']
    arrtab.rows = ['start', 'mid', 'end']
    print(arrtab.headers())
    print(arrtab.col(0))
    print(arrtab.col('Step'))
    print(arrtab.row(1))  # index access
    print(arrtab.row('mid'))  # name access
    print(arrtab.rowcol('end', 'Temp'))
    print(arrtab.rowcol(2, 1))
    arrtab.console(delim=',')

    print()
    rowtab = ListListTable([[0, 300], [1, 475], [2, 600]])
    print(rowtab.shape())
    rowtab.keywords = ['Step', 'Temp']
    rowtab.rows = ['start', 'mid', 'end']
    print(rowtab.headers())
    print(rowtab.col(0))
    print(rowtab.col('Step'))
    print(rowtab.row(1))  # index access
    print(rowtab.row('mid'))  # name access
    print(rowtab.rowcol('end', 'Temp'))
    print(rowtab.rowcol(2, 1))

    rowtab.console(delim='|', spacer='-' * 9)

    rowtab.col_fun('sum')
    rowtab.col_fun('avg')
    rowtab.col_fun('avg0')  # avg ignoring 0's

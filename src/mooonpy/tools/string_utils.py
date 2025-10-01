# -*- coding: utf-8 -*-
import re
from typing import Union
import numpy as np

def is_float(string: str) -> bool:
    """
    regex to check if a string is a float
    """
    float_re = re.compile(r'^-?\d+(\.\d+)?([eE][-+]?\d+)?$')
    return bool(float_re.match(string))

def string2digit(string: str) -> Union[int, str, float]:
    """
    Converts a string to float or int
    :param string: String to be converted
    :type string: str
    :return: int or float or sring

    :Example:
        >>> from mooonpy.tools import string2digit
        >>> print(string2digit('5'))
        5
        >>> print(string2digit('5.1'))
        5.1
        >>> print(string2digit('5a'))
        '5a'
    """
    if string.isnumeric():
        return int(string)
    elif is_float(string):
        return float(string)
    else:
        return string

def list2digit(string_list):
    return [string2digit(string) for string in string_list]

def _col_convert(column,skip_int=False):
    try:
        column = np.array(column, float)  ## convert from string. This is about half the runtime
    except:
        raise Exception(f'Column not convertable to floats: {column}')
    if skip_int: return column # cannot convert to int, throws warning if the next line executes with a nan
    col_int = column.astype(int)  # convert to int array. exact up to 9 quadrillion
    if np.all(column == col_int): # can't find a faster way to do this, it gets all truth values first, but it really should have a breakout for first float found, but the C functions are fast
        return col_int
    else:
        return column
    # for f_col,i_col in zip(column,col_int): # not faster, probably loop setup
    #     if f_col != i_col:
    #         return column
    # return col_int

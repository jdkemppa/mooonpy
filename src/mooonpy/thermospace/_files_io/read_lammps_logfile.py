from mooonpy import Path, Thermospace
from mooonpy.tools.string_utils import _col_convert
import warnings


def read_lammps_logfile(file: Path | str, silence_error_line: bool = False, out: Thermospace = None) -> Thermospace:
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
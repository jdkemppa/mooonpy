def readcsv_basic(file: path):
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
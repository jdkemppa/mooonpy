# -*- coding: utf-8 -*-
from bz2 import open as bz2_open
from glob import glob
from gzip import open as gzip_open
from lzma import open as lzma_open
from os import sep
from os.path import normpath, join, exists, abspath, basename, dirname, splitext, getmtime
from typing import List, Optional, Union


class Path(str):
    """
    *As computational scientists, half our jobs is file management and manipulation,
    the Path class contains several aliases for the os.path and glob.glob modules
    to make processing data easier. All mooonpy functions internally use this class
    for inputs of files or folders. Relevant strings are converted to path on entering functions*

    Examples
    --------
    A copy of the code used in these examples is avalible in root\\mooonpy\\examples\\tools\\path_utils\\example_Path.py

    **Basic Path Operations**
        >>> project_path = Path('Project/Data/Analysis')
        >>> filename = Path('results.txt')
        >>> full_path = project_path / filename
        >>> print(full_path)
        Project\\Data\\Analysis\\results.txt
        >>> print(abs(full_path))
        root\\mooonpy\\examples\\tools\\path_utils\\Project\\Data\\Analysis\\results.txt

    **Path Parsing**
        >>> sample_path = Path('experiments/run_001/data.csv.gz')
        >>> print(sample_path.dir())
        experiments\\run_001
        >>> print(sample_path.basename())
        data.csv.gz
        >>> print(sample_path.root())
        data.csv
        >>> print(sample_path.ext())
        .gz

    **Extension Manipulation**
        >>> data_file = Path('analysis/results.txt')
        >>> print(data_file.new_ext('.json'))
        analysis\\results.json
        >>> print(data_file.new_ext('.txt.gz'))
        analysis\\results.txt.gz

    **File Existence**
        >>> current_file = Path(__file__)
        >>> fake_file = Path('nonexistent.txt')
        >>> print(bool(current_file))
        True
        >>> print(bool(fake_file))
        False

    **Wildcard Matching**
        >>> txt_pattern = Path('temp_dir/*.txt')
        >>> print(txt_pattern.matches())
        ['test1.txt', 'test2.txt']
        >>> for file in Path('temp_dir/*'):
        ...     print(file.basename())
        data.csv
        readme.md
        test1.txt
        test2.txt

    **Recent File Finding**
        >>> pattern = Path('temp_dir/*.txt')
        >>> print(pattern.recent())
        newest_file.txt
        >>> print(pattern.recent(oldest=True))
        old_file.txt

    **Smart File Opening**
        >>> mypath = Path('data.txt')
        >>> with mypath.open('w') as f:
        ...     f.write('Hello World')
        # Creates regular file
        >>> compressed_path = Path('data.txt.gz')
        >>> # compressed_path.open() would use gzip automatically
        # Would automatically handle gzip compression
        ** Absolute Path Conversion **
        >>> rel_path = Path('data/file.txt')
        >>> print(abs(rel_path))
        root\\mooonpy\\examples\\tools\\path_utils\\data\\file.txt

    .. TODO::
        __truediv__ __bool__ __abs__ and __iter__ docstrings in config?
    """
    search_prefixes = None  # Optional[List['Path']] — full paths up to and including common
    search_common = None  # Optional[str]          — common dir name for splitting abs paths

    # Column-wise natsort direction for __iter__ over wildcard matches.
    #   True  — sort right-to-left: rightmost wildcard is the tiebreaker,
    #           leftmost is the primary key (default).
    #   False — sort left-to-right: leftmost wildcard is the tiebreaker,
    #           rightmost is the primary key.
    #   None  — disable natsort entirely; iterate in matches() order.
    natsort_right_to_left = True
    #   If True, drop matches whose value in the highest-priority wildcard
    #   column is non-numeric (e.g. 'last', 'crashed'). Useful when such
    #   matches are duplicates of a numbered run and should be excluded from
    #   analysis. Default False — keep them, sorted after numeric values.
    ignore_nonnumeric_priority_wild = False

    def __fspath__(self) -> str:
        return str(self)  # Mostly fixes type hints

    def __new__(cls, string: Union[str, 'Path']) -> 'Path':  ## this is before init somehow
        return super().__new__(cls, normpath(string))  # Typing is confused here

    def __truediv__(self, other: Union[str, 'Path']):
        """
        Join paths with subdirectory delimiter (ie /).

        Alias for os.path.join.

        :param other: Path to join on right
        :type other: str or Path
        :param self: Path to join on left
        :return: Joined Path
        :rtype: Path

        :Example:
            >>> from mooonpy.tools import Path
            >>> MyDir = Path('Project/Monomers')
            >>> MyFile = Path('DETDA.mol')
            >>> print(MyDir / MyFile)
            'Project\\Monomers\\DETDA.mol'
        """
        # return Path(join(self, other)) # does not work in Linux or Mac
        return Path(join(str(self), str(other)))  # fixes

    def __sub__(self, template: Union[str, 'Path']):
        """
        Returns value of glob wildcard (*) characters from a filename-template pair using '-'
        
        :param other: template with glob wildcards
        :type template: str or Path with wildcards
        :param self: Path or str without wildcards
        :return: list of wildcards
        :rtype: list

        :Example:
          >>> from mooonpy.tools import Path
          >>> myfilename = Path('filename_step10_temp300.data')
          >>> mytemplate = Path('filename_step*_temp*.data')
          >>> print(myfilename - mytemplate)
          ['10', '300']
        """
        # import built-in matching tools
        import fnmatch, re

        def match(self, template):
            # convert paths to strings for easier treatment
            real_file = str(self).replace('/', '\\')  # fix dash for linux compatability
            template_file = str(template).replace('/', '\\')
            # check if template matches the tested file
            if fnmatch.fnmatch(real_file, template_file):
                # build a regrex expression, then pull a re match out
                to_match = fnmatch.translate(template_file)  # convert to regular expression
                to_match = to_match.replace('\\Z(?ms)', r'\Z')
                to_match = to_match.replace('.*', '(.*?)')  # allow un-greedy capturable regrex
                compiled = re.compile(to_match, re.S | re.IGNORECASE)  # compile regrex string into object
                matched = compiled.match(real_file)
                return matched
            else:
                # return nothing if no match was found
                return None

                # try matching path to template

        matched = match(self, template)
        if matched == None:
            longer_matched = match(self, Path('*') / template)
            if longer_matched == None:
                return None
            else:
                return list(longer_matched.groups())  # try again, prepending a */
        else:
            return list(matched.groups())

    def __bool__(self) -> bool:
        """
        Check if the path points to a file or directory.

        Alias for os.path.exists.

        :return: Boolean True/False
        :rtype: bool

        :Example:
            >>> from mooonpy.tools import Path
            >>> MyFile1 = Path('DETDA.mol')
            >>> print(is MyFile1)
            True
            >>> MyFile2 = Path('doesnotexist.mol')
            >>> print(is MyFile2)
            False
        """
        return exists(str(self))

    def __abs__(self) -> 'Path':
        """
        Absolute path to file or directory.

        Alias for os.path.abspath.

        :return: Absolute path to file or directory
        :rtype: Path

        :Example:
            >>> from mooonpy import Path
            >>> MyFile = Path('DETDA.mol')
            >>> print(abs(MyFile))
            'C:\\Users\\You\\Desktop\\DETDA.mol'
        """
        return Path(abspath(self))

    def __iter__(self) -> iter:
        """
        Iterates through matching paths with a * (asterisk) wildcard character.

        Matches are sorted column-by-column over the captured wildcard values
        (as extracted by :meth:`__sub__`). The traversal direction is
        controlled by :attr:`natsort_right_to_left`:

        * ``True`` *(default)* — right-to-left: the rightmost wildcard is the
          primary sort key, the leftmost is the final tiebreaker. This
          matches the common case of ``run_*/step_*.data`` where the
          trailing wildcard (``step``) is the dominant integer counter.
        * ``False`` — left-to-right: the leftmost wildcard is the primary
          sort key and the rightmost is the final tiebreaker.
        * ``None`` — disable natsort; iterate in the order produced by
          :meth:`matches` (lexicographic, from :func:`glob.glob`).

        Within each column, values are sorted with all numeric values
        (signed ints and floats) first in numeric order, followed by
        non-numeric strings in lexicographic order. This lets mixed
        columns like step values of ``1, 2, 10, 'last', 'crashed'``
        sort cleanly as ``1, 2, 10, 'crashed', 'last'`` without
        crashing or falling back to pure lex order. Blank captures sort
        as the last float (via ``inf``) and the first non-numeric string
        (empty string is lex-minimum).

        If :attr:`ignore_nonnumeric_priority_wild` is ``True``, matches whose
        value in the primary column (rightmost when ``natsort_right_to_left``
        is ``True``, leftmost when ``False``) is non-numeric are dropped
        entirely. This is useful when names like ``step_last.data`` are
        duplicates of a numbered run and should be excluded from analysis.

        :return: iter object of List of matching Paths
        :rtype: iter

        .. note:: This overrides string iteration through characters; convert
            back to string before passing into a function if this causes issues.
        .. note:: The natsort method used may not behave well if wildcards intended
            as strings may be interpreted as floats, where file_nan.data
            from file_*.data will be sorted as a float not string
        :Example:
            >>> from mooonpy import Path
            >>> for f in Path('temp_dir/output_*.data'):
            ...     print(f.basename())
            output_1.data
            output_2.data
            output_10.data
            output_crashed.data
            output_last.data

            >>> Path.ignore_nonnumeric_priority_wild = True
            >>> for f in Path('temp_dir/output_*.data'):
            ...     print(f.basename())
            output_1.data
            output_2.data
            output_10.data

            >>> # Multiple wildcards, right-to-left (default):
            >>> # primary key is step (rightmost), run breaks ties.
            >>> for f in Path('run_*/step_*.data'):
            ...     print(f)
        """
        found = self.matches()
        if (self.natsort_right_to_left is None
                or '*' not in str(self)
                or len(found) < 2):
            return iter(found)

        # Pair each match with its captured wildcard tuple.
        rows = []
        for m in found:
            groups = m - self # call __sub__ operator
            if not groups:
                from warnings import warn
                warn('Path __sub__ operator failed in __iter__ operator call, returned unsorted matches')
                return iter(found)
            rows.append((tuple(groups), m))

        ncols = len(rows[0][0])
        if not all(len(r[0]) == ncols for r in rows):
            from warnings import warn
            warn('Path found uneven wildcards in __iter__ operator call, returned unsorted matches')
            return iter(found)  # ragged captures; bail safely

        # Column priority order: right-to-left means we start at the rightmost
        # column (it's the primary key), then walk leftward for tiebreakers.
        col_order = (list(range(ncols - 1, -1, -1))
                     if self.natsort_right_to_left
                     else list(range(ncols)))

        def as_number(v):
            if v == '':
                return float('inf')
            try:
                num = float(v)
            except ValueError:
                return None
            return num

        def col_key(value: str):
            # (0, num) for numeric values, (1, str) for non-numeric.
            # Tuple ordering puts all numeric before all non-numeric, and
            # within each group the second element orders them naturally.
            num = as_number(value)
            return (0, num) if num is not None else (1, value)

        # Optionally drop rows whose primary-column value is non-numeric.
        if self.ignore_nonnumeric_priority_wild:
            primary = col_order[0]
            rows = [r for r in rows if as_number(r[0][primary]) is not None]
            if len(rows) < 2:
                return iter([m for _, m in rows])

        # Recursively sort: fix higher-priority columns first, then sort each
        # group by the next column with numeric-before-string ordering.
        def sort_rows(subset, depth):
            if depth >= len(col_order) or len(subset) < 2:
                return subset
            col = col_order[depth]
            subset = sorted(subset, key=lambda r: col_key(r[0][col]))
            # Group by this column's value (post-sort) and recurse into each.
            result = []
            i = 0
            while i < len(subset):
                j = i + 1
                key_i = subset[i][0][col]
                while j < len(subset) and subset[j][0][col] == key_i:
                    j += 1
                result.extend(sort_rows(subset[i:j], depth + 1))
                i = j
            return result

        return iter([m for _, m in sort_rows(rows, 0)])

    def basename(self) -> 'Path':
        """
        Split Path to filename and extention.

        Alias for os.path.basename

        :return: Path of file
        :rtype: Path

        :Example:
            >>> from mooonpy import Path
            >>> MyPath = Path('Project/Monomers/DETDA.mol')
            >>> print(MyPath.basename())
            'DETDA.mol'
        """
        return Path(basename(self))

    def dir(self) -> 'Path':
        """
        Split Path to directory.

        Alias for os.path.dirname.

        :return: Path to directory
        :rtype: Path

        :Example:
            >>> from mooonpy import Path
            >>> MyPath = Path('Project/Monomers/DETDA.mol')
            >>> print(MyPath.dir())
            'Project\\Monomers'
        """
        return Path(dirname(self))

    def ext(self) -> 'Path':
        """
        Split Path to just extention.

        Alias for os.path.basename and splitext.

        :return: extention as Path
        :rtype: Path

        :Example:
            >>> from mooonpy import Path
            >>> MyPath = Path('Project/Monomers/DETDA.mol')
            >>> print(MyPath.ext())
            '.mol'
        """
        return Path(splitext(self.basename())[1])

    def matches(self, whitelist_ext=None, blacklist_ext=None) -> List['Path']:
        """
        Finds matching paths with a * (asterisk) wildcard character.

        :return: List of matching Paths
        :rtype: List[Path]

        :Example:
            >>> from mooonpy import Path
            >>> MyWildcard = Path('*.mol')
            >>> print(Path.matches(MyWildcard))
            [Path('DETDA.mol'), Path('DEGBF.mol')]
        """
        matches = [Path(file) for file in glob(self)]
        if blacklist_ext:
            matches = [match for match in matches if match.ext() not in blacklist_ext]
        elif whitelist_ext:
            matches = [match for match in matches if match.ext() in whitelist_ext]

        return matches

    def new_ext(self, ext: Union[str, 'Path']) -> 'Path':
        """
        Replace extension on a Path with a new extension.

        :param ext: new extension including delimeter.

        :type ext: str or Path
        :return: replaced Path
        :rtype: Path

        :Example:
            >>> from mooonpy import Path
            >>> MyPath = Path('Project/Monomers/DETDA.mol')
            >>> print(MyPath.new_ext('.data'))
            'Project/Monomers/DETDA.data'
        """
        return Path(splitext(self)[0] + ext)

    def open(self, mode='r', encoding='utf-8'):
        """
        Open path with smart_open

        :param mode: Open mode, usually 'r' or 'a'
        :type mode: str
        :param encoding: File encoding
        :type encoding: str
        :return: opened file as object
        :rtype: File Object

        :Example:
            >>> from mooonpy import Path
            >>> MyPath = Path('Project/Monomers/DETDA.mol')
            >>> MyFileObj = MyPath.open(mode='r')

        """
        return smart_open(self, mode, encoding)

    def recent(self, oldest: bool = False) -> Optional['Path']:
        """
        Find wildcard matches and return the Path of the most recently modified file.

        :param oldest: Reverses direction and finds least recently modified file.
        :type oldest: bool
        :return: Path of most recently modified file
        :rtype: Path

        :Example:
            >>> from mooonpy import Path
            >>> MyWildcard = Path('Template_*.lmpmol')
            >>> print(Path.recent())
            'Template_1_v10_final_realthistime.lmpmol'
            >>> print(Path.recent(oldest=True))
            'Template_1.lmpmol'
        """
        times = {}
        if self.search_prefixes:
            for file in self.locate_all():
                times[getmtime(file)] = file
        else:
            for file in self:
                times[getmtime(file)] = file
        if times:
            sorted_time = sorted(list(times.keys()))
            if oldest:
                return times[sorted_time[0]]
            else:
                return times[sorted_time[-1]]
        else:
            return None

    def root(self) -> 'Path':
        """
        Split Path to filename with no extention.

        Alias for os.path.basename and splitext.

        :return: Path of filename
        :rtype: Path

        :Example:
            >>> from mooonpy import Path
            >>> MyPath = Path('Project/Monomers/DETDA.mol')
            >>> print(MyPath.root())
            'DETDA'
        """
        return Path(splitext(self.basename())[0])

    def format(self, *args, **kwargs):
        return Path(str(self).format(*args, **kwargs))

    @classmethod
    def find_prefix(cls, common=None, path=None, add=True) -> Optional['Path']:
        """
        Extract the prefix of *path* up to and including *common*, and
        optionally append it to :attr:`search_prefixes`.

        When *path* is omitted the caller's ``__file__`` is used automatically,
        so a script located inside the common directory tree only needs to supply
        the common directory name.  When *common* is omitted the already-stored
        :attr:`search_common` is reused, allowing multiple calls for different
        drives without repeating the directory name.

        :param common: Shared directory name, e.g. ``'research'``.  Omit to
                       reuse :attr:`search_common` (must have been set earlier).
        :type common:  str or None
        :param path:   Path containing *common* as a component.  Omit to use
                       the calling script's ``__file__`` automatically.
        :type path:    str, Path, or None
        :param add:    Append the extracted prefix to :attr:`search_prefixes`
                       (default ``True``).  Duplicates are skipped.
        :type add:     bool
        :return:       Extracted prefix ``Path``, or ``None`` if *common* was
                       not found in *path*.
        :rtype:        Path or None

        :Example:
            >>> # In a script at C:/research/sims/run_001/analysis.py:
            >>> Path.find_prefix('research')           # auto-detects C:\\research
            Path('C:\\\\research')
            >>> Path.find_prefix(path='D:/backups/research/sims/run_001/analysis.py')
            Path('D:\\\\backups\\\\research')          # reuses search_common='research'
        """
        import inspect
        if path is None:
            path = Path(inspect.stack()[1].filename)
        else:
            path = Path(path)
        if common is None:
            common = cls.search_common  # reuse previously set value
        if common is None:
            raise RuntimeError(
                "Path.search_common is not set — pass 'common' explicitly or call "
                "Path.find_prefix() with a common directory name first."
            )
        s = normpath(str(path))
        before, found, _ = s.partition(sep + str(common) + sep)
        if not found:
            if s.endswith(sep + str(common)):  # common is the last component
                prefix = Path(s)
            else:
                return None
        else:
            prefix = Path(before + sep + str(common))
        if cls.search_common is None:
            cls.search_common = str(common)
        if add:
            if cls.search_prefixes is None:
                cls.search_prefixes = []
            if prefix not in cls.search_prefixes:
                cls.search_prefixes.append(prefix)
        return prefix

    def _abspath_from_common(self, prefix) -> 'Path':
        """Return *prefix* joined with the portion of self that follows search_common.

        Resolves *self* to an absolute path via CWD first, so a relative path
        passed from a script inside the common hierarchy automatically carries
        its subdirectory context to the new prefix.

        If common is not found in the resolved path (e.g. a plain relative tail
        like ``'sims/run_001/out.log'`` used directly), *self* is treated as
        already being the tail and is appended to *prefix* unchanged.

        :param prefix: Target prefix to prepend (str or Path).
        :return:       Absolute rewritten Path under *prefix*.
        :rtype:        Path
        """
        if Path.search_common is None:
            raise RuntimeError(
                "Path.search_common is not set — cannot strip the drive prefix.  "
                "Call Path.find_prefix() first."
            )
        # Always resolve to absolute so relative paths inherit CWD context.
        s = abs(self)
        needle = sep + str(Path.search_common) + sep
        needle_end = sep + str(Path.search_common)
        _, found, tail = s.partition(needle)
        if found:
            return Path(prefix) / tail
        if s.endswith(needle_end):  # common is the last component; no tail
            return Path(prefix)
        # common not found — treat self as an already-relative tail from common
        return Path(prefix) / self

    def swap_prefix(self, target) -> Optional['Path']:
        """
        Return a new path with this path's prefix replaced by *target*.

        Strips everything up to and including :attr:`search_common` from *self*
        (via :meth:`_relpath_from_common`), then prepends the chosen prefix.

        :param target: Replacement prefix — either an ``int`` index into
                       :attr:`search_prefixes`, or a ``str``/``Path`` value.
        :type target:  int, str, or Path
        :return:       Rewritten ``Path``, or ``None`` if *target* is an
                       out-of-range integer index.
        :rtype:        Path or None

        :Example:
            >>> Path.search_prefixes = [Path('C:/research'), Path('D:/backups/research')]
            >>> Path.search_common   = 'research'
            >>> p = Path('D:/backups/research/sims/run_001/output.log')
            >>> p.swap_prefix(0)
            Path('C:\\\\research\\\\sims\\\\run_001\\\\output.log')
            >>> p.swap_prefix('D:/backups/research')
            Path('D:\\\\backups\\\\research\\\\sims\\\\run_001\\\\output.log')
        """
        if isinstance(target, int):
            if Path.search_prefixes is None or target >= len(Path.search_prefixes):
                return None
            prefix = Path.search_prefixes[target]
        else:
            prefix = Path(target)
        return self._abspath_from_common(prefix)

    def locate(self, recent=None) -> Optional['Path']:
        """
        Locate this path across :attr:`search_prefixes`.

        Behaviour depends on *recent* and whether the path contains a ``*``
        wildcard:

        * ``recent=None`` *(default)* — **prefix-order priority**

          - *No wildcard*: return the first prefix for which the exact file
            exists.
          - *Wildcard*: return the glob **pattern** (``*`` kept in place) for
            the first prefix that has any matches.  Pass that result to
            :meth:`matches` or iterate over it to expand the files.

        * ``recent=True`` — return the **most recently modified** file across
          all prefixes (wildcards expanded via :meth:`locate_all`).
        * ``recent=False`` — return the **oldest** file across all prefixes.

        Naming note: ``locate`` avoids shadowing the built-in
        :meth:`str.find` method.

        :param recent: ``None`` for prefix-order (default), ``True`` for
                       newest, ``False`` for oldest.
        :return:       Matched ``Path``, or ``None``.
        :rtype:        Path or None

        :Example:
            >>> Path.search_prefixes = [Path('C:/research'), Path('D:/backups/research')]
            >>> Path.search_common   = 'research'
            >>> Path('sims/run_001/output.log').locate()
            Path('C:\\\\research\\\\sims\\\\run_001\\\\output.log')
            >>> Path('sims/*/output.log').locate()          # returns the pattern
            Path('C:\\\\research\\\\sims\\\\*\\\\output.log')
            >>> Path('sims/*/output.log').locate(recent=True)
            Path('D:\\\\backups\\\\research\\\\sims\\\\run_002\\\\output.log')
            >>> Path('sims/*/output.log').locate(recent=False)
            Path('C:\\\\research\\\\sims\\\\run_001\\\\output.log')
        """
        if Path.search_prefixes is None:
            raise RuntimeError(
                "Path.search_prefixes is not configured — call Path.find_prefix() first."
            )
        if not Path.search_prefixes:  # configured but empty list → no matches possible
            return None
        if recent is not None:
            all_matches = self.locate_all()
            if not all_matches:
                return None
            key = lambda p: getmtime(str(p))
            return max(all_matches, key=key) if recent else min(all_matches, key=key)
        has_wildcard = '*' in str(self)
        for i in range(len(Path.search_prefixes)):
            candidate = self.swap_prefix(i)
            if candidate is None:
                continue
            if has_wildcard:
                if candidate.matches():  # any glob hits on this drive?
                    return candidate  # return the pattern (still has *)
            else:
                if candidate:  # exact file exists?
                    return candidate
        return None

    def locate_all(self, whitelist_ext=None, blacklist_ext=None) -> List['Path']:
        """
        Return all glob matches of this path across every prefix in
        :attr:`search_prefixes`.

        Wildcards are expanded on each prefix independently, so
        ``Path('sims/*/output.log').locate_all()`` collects every matching file
        across all configured drives.

        :param whitelist_ext: If given, only include paths with these extensions.
        :param blacklist_ext: If given, exclude paths with these extensions.
        :return:  All matching ``Path`` objects across all prefixes.
        :rtype:   List[Path]

        :Example:
            >>> Path.search_prefixes = [Path('C:/research'), Path('D:/backups/research')]
            >>> Path.search_common   = 'research'
            >>> Path('sims/*/output.log').locate_all()
            [Path('C:\\\\research\\\\sims\\\\run_001\\\\output.log'),
             Path('D:\\\\backups\\\\research\\\\sims\\\\run_001\\\\output.log'),
             Path('D:\\\\backups\\\\research\\\\sims\\\\run_002\\\\output.log')]
        """
        if Path.search_prefixes is None:
            raise RuntimeError(
                "Path.search_prefixes is not configured — call Path.find_prefix() first."
            )
        if not Path.search_prefixes:  # configured but empty list
            return []
        results = []
        for i in range(len(Path.search_prefixes)):
            pattern = self.swap_prefix(i)
            if pattern is not None:
                results.extend(pattern.matches(whitelist_ext=whitelist_ext,
                                               blacklist_ext=blacklist_ext))
        return results


# End of Path

# %% Misc file tools
def smart_open(filename, mode='r', encoding='utf-8'):
    """
    Open file with appropriate decompression based on extension

    **Supported extensions: Use substring in filename**
        - .gz: Uses gzip module
        - .bz2: Uses bzip2 module
        - .xz: Uses lzma module
        - .lzma: Uses lzma module
        - Other extensions use the builtin open function


    :param filename: Path to file
    :type filename: Path or str
    :param mode: Open mode, usually 'r', 'w' or 'a'
    :type mode: str
    :param encoding: File encoding
    :type encoding: str

    :return: opened file as object
    :rtype: File Object
    :Example:
        >>> from mooonpy.tools.file_utils import smart_open
        >>> MyFileObj = smart_open('Project/Monomers/DETDA.data.gz')
    """
    try:
        if '.gz' in filename:
            return gzip_open(str(filename), mode + 't', encoding=encoding)
        elif '.bz2' in filename:
            return bz2_open(str(filename), mode + 't', encoding=encoding)
        elif '.xz' in filename or '.lzma' in filename:
            return lzma_open(str(filename), mode + 't', encoding=encoding)
    except:

        pass  # compressed filename did not work
    return open(str(filename), mode, encoding=encoding)  # try regular read

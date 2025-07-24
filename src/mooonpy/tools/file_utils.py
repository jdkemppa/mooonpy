# -*- coding: utf-8 -*-
import bz2
import glob
import gzip
import lzma
import os
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
    def __fspath__(self) -> str:
        return str(self)  # Mostly fixes type hints

    def __new__(cls, string: Union[str, 'Path']) -> 'Path':  ## this is before init somehow
        return super().__new__(cls, os.path.normpath(string))  # Typing is confused here

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
        # return Path(os.path.join(self, other)) # does not work in Linux or Mac
        return Path(os.path.join(str(self), str(other))) # fixes

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
        return os.path.exists(self)

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
        return Path(os.path.abspath(self))

    def __iter__(self) -> iter:
        """
        Iterates through matching paths with a * (asterisk) wildcard character.

        :return: iter object of List of matching Paths

        :rtype: iter

        .. note:: This overrides string iteration through characters, convert back to string
        before passing into a function if this causes issues.

        :Example:
            >>> from mooonpy import Path
            >>> MyWildcard = Path('*.mol')
            >>> for MyMatch in MyWildcard:
            >>>     print(MyMatch)
            'DETDA.mol'
            'DEGBF.mol'
        """
        return iter(self.matches())

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
        return Path(os.path.basename(self))

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
        return Path(os.path.dirname(self))

    def ext(self) -> 'Path':
        """
        Split Path to just extention.

        Alias for os.path.basename and os.path.splitext.

        :return: extention as Path
        :rtype: Path

        :Example:
            >>> from mooonpy import Path
            >>> MyPath = Path('Project/Monomers/DETDA.mol')
            >>> print(MyPath.ext())
            '.mol'
        """
        return Path(os.path.splitext(self.basename())[1])

    def matches(self,whitelist_ext=None, blacklist_ext=None) -> List['Path']:
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
        matches = [Path(file) for file in glob.glob(self)]
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
        return Path(os.path.splitext(self)[0] + ext)

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
        for file in self:
            times[os.path.getmtime(file)] = file
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

        Alias for os.path.basename and os.path.splitext.

        :return: Path of filename
        :rtype: Path

        :Example:
            >>> from mooonpy import Path
            >>> MyPath = Path('Project/Monomers/DETDA.mol')
            >>> print(MyPath.root())
            'DETDA'
        """
        return Path(os.path.splitext(self.basename())[0])

# End of Path

#%% Misc file tools
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
            return gzip.open(str(filename), mode + 't', encoding=encoding)
        elif '.bz2' in filename:
            return bz2.open(str(filename), mode + 't', encoding=encoding)
        elif '.xz' in filename or '.lzma' in filename:
            return lzma.open(str(filename), mode + 't', encoding=encoding)
    except:

        pass  # compressed filename did not work
    return open(str(filename), mode, encoding=encoding)  # try regular read


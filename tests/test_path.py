# -*- coding: utf-8 -*-

import pytest
from mooonpy.tools.file_utils import Path, smart_open

import os
import tempfile
import gzip
import bz2
import lzma
import shutil
from pathlib import Path as StdPath

# Run pytest tests
print("=== NOT Running Tests ===\n")
print("To run the pytest tests, use:")
print("   Command: pytest test_path.py -v")
print("   Command: pytest test_path.py::TestPath::test_path_creation -v  # Run specific test")
print("   Command: pytest test_path.py::TestSmartOpen -v  # Run specific test class")
print()

class TestPath:
    """Pytest tests for Path class"""

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory with test files"""
        temp_dir = Path(tempfile.mkdtemp())
        # temp_dir = Path(__file__).dir() / "temp_dir"
        # os.makedirs(temp_dir, exist_ok=True)
        test_files = []

        # Create test files
        for i, name in enumerate(["file1.txt", "file2.txt", "data.csv", "script.py"]):
            filepath = os.path.join(temp_dir, name)
            with open(filepath, 'w') as f:
                f.write(f"Content {i}")
            test_files.append(filepath)

            # Add different modification times
            import time
            time.sleep(0.01)  # Small delay to ensure different mtimes

        yield temp_dir, test_files

        # Cleanup
        shutil.rmtree(temp_dir)

    def test_path_creation(self):
        """Test Path object creation and normalization"""
        p1 = Path("folder/file.txt")
        p2 = Path("folder//file.txt")  # Double slash
        p3 = Path("folder\\file.txt")  # Backslash

        # All should be normalized
        assert isinstance(p1, Path)
        assert isinstance(p1, str)
        assert str(p1) == os.path.normpath("folder/file.txt")

    def test_path_division(self):
        """Test path joining with / operator"""
        base = Path("project")
        sub = Path("data")
        filename = "file.txt"

        result = base / sub / filename
        expected = Path(os.path.join("project", "data", "file.txt"))

        assert result == expected
        assert isinstance(result, Path)

    def test_path_existence(self, temp_dir):
        """Test file existence checking"""
        temp_dir_path, test_files = temp_dir
        existing_file = Path(test_files[0])
        non_existing = Path("nonexistent_file.txt")

        assert bool(existing_file) == True
        assert bool(non_existing) == False

    def test_absolute_path(self):
        """Test absolute path conversion"""
        relative = Path("test.txt")
        absolute = abs(relative)

        assert os.path.isabs(absolute)
        assert isinstance(absolute, Path)

    def test_path_parsing(self):
        """Test path component extraction"""
        test_path = Path("folder/subfolder/file.txt")

        assert test_path.basename() == "file.txt"
        assert test_path.dir() == Path("folder/subfolder")
        assert test_path.root() == "file"
        assert test_path.ext() == ".txt"

    def test_extension_replacement(self):
        """Test extension replacement"""
        original = Path("data/file.txt")
        new_path = original.new_ext(".csv")

        assert new_path == Path("data/file.csv")
        assert isinstance(new_path, Path)

    def test_wildcard_matching(self, temp_dir):
        """Test wildcard pattern matching"""
        temp_dir_path, test_files = temp_dir
        pattern = Path(os.path.join(temp_dir_path, "*.txt"))
        matches = pattern.matches()

        assert len(matches) == 2  # file1.txt and file2.txt
        assert all(isinstance(m, Path) for m in matches)
        assert all(m.ext() == ".txt" for m in matches)

    def test_iterator(self, temp_dir):
        """Test Path iteration over wildcard matches"""
        temp_dir_path, test_files = temp_dir
        pattern = Path(os.path.join(temp_dir_path, "*"))
        files = list(pattern)

        assert len(files) == 4  # All test files
        assert all(isinstance(f, Path) for f in files)

    def test_recent_file(self, temp_dir):
        """Test finding most/least recently modified files"""
        temp_dir_path, test_files = temp_dir
        pattern = Path(os.path.join(temp_dir_path, "*.txt"))

        most_recent = pattern.recent()
        oldest = pattern.recent(oldest=True)

        assert isinstance(most_recent, str)
        assert isinstance(oldest, str)
        assert most_recent != oldest

    def test_recent_no_matches(self, temp_dir):
        """Test recent() with no matching files"""
        temp_dir_path, test_files = temp_dir
        pattern = Path(os.path.join(temp_dir_path, "*.nonexistent"))
        result = pattern.recent()

        assert result is None


class TestSmartOpen:
    """Pytest tests for smart_open function"""

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for tests"""
        temp_dir = Path(tempfile.mkdtemp())
        # temp_dir = Path(__file__).dir() / "temp_dir"
        # os.makedirs(temp_dir, exist_ok=True)
        yield temp_dir
        shutil.rmtree(temp_dir)

    @pytest.fixture
    def test_content(self):
        """Test content for file operations"""
        return "Hello, World!\nTest content here."

    def test_regular_file(self, temp_dir, test_content):
        """Test opening regular files"""
        filepath = temp_dir / "test.txt"

        # Write file
        with smart_open(filepath, 'w') as f:
            f.write(test_content)

        # Read file
        with smart_open(filepath, 'r') as f:
            content = f.read()

        assert content == test_content

    def test_gzip_file(self, temp_dir, test_content):
        """Test opening gzip files"""
        filepath = temp_dir / "test.txt.gz"

        # Write compressed file
        with smart_open(filepath, 'w') as f:
            f.write(test_content)

        # Read compressed file
        with smart_open(filepath, 'r') as f:
            content = f.read()

        assert content == test_content

    def test_bzip2_file(self, temp_dir, test_content):
        """Test opening bzip2 files"""
        filepath = temp_dir / "test.txt.bz2"

        # Write compressed file
        with smart_open(filepath, 'w') as f:
            f.write(test_content)

        # Read compressed file
        with smart_open(filepath, 'r') as f:
            content = f.read()

        assert content == test_content

    def test_lzma_file(self, temp_dir, test_content):
        """Test opening LZMA files"""
        for ext in ['.xz', '.lzma']:
            filepath = temp_dir / f"test.txt{ext}"

            # Write compressed file
            with smart_open(filepath, 'w') as f:
                f.write(test_content)

            # Read compressed file
            with smart_open(filepath, 'r') as f:
                content = f.read()

            assert content == test_content

    # def test_fallback_to_regular_open(self, temp_dir, test_content):
    #     """Test fallback when compression fails"""
    #     filepath = temp_dir / "test.txt.gz"
    #
    #     # Create a file that looks compressed but isn't
    #     with open(filepath, 'w') as f:
    #         f.write(test_content)
    #
    #     # Should still be able to read it.
    #     # No, it returns the file handle but fails on read. Error message is complete enough
    #     with smart_open(filepath, 'r') as f:
    #         content = f.read()
    #
    #     assert content == test_content


class TestMultiDrivePath:
    """Tests for multi-drive/multi-prefix: find_prefix, swap_prefix, locate, locate_all."""

    @pytest.fixture(autouse=True)
    def reset_class_vars(self):
        """Reset Path class variables before and after every test."""
        Path.search_prefixes = None
        Path.search_common = None
        yield
        Path.search_prefixes = None
        Path.search_common = None

    @pytest.fixture
    def two_drives(self):
        """Two temp dirs simulating two drive prefixes sharing a common subdir."""
        drive1 = Path(tempfile.mkdtemp())
        drive2 = Path(tempfile.mkdtemp())
        common = 'research'

        prefix1 = drive1 / common
        prefix2 = drive2 / common
        os.makedirs(prefix1, exist_ok=True)
        os.makedirs(prefix2, exist_ok=True)

        yield prefix1, prefix2

        shutil.rmtree(drive1)
        shutil.rmtree(drive2)

    def _make_file(self, prefix, relpath):
        """Create a file at prefix/relpath, making intermediate dirs."""
        full = Path(prefix) / relpath
        os.makedirs(str(Path(full).dir()), exist_ok=True)
        with open(full, 'w') as f:
            f.write('test')
        return full

    # ------------------------------------------------------------------
    # locate() — exact paths
    # ------------------------------------------------------------------

    def test_locate_first_prefix(self, two_drives):
        """locate() returns file when it exists only on the first prefix."""
        prefix1, prefix2 = two_drives
        relpath = 'sims/run_001/output.log'
        expected = self._make_file(prefix1, relpath)

        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        result = Path(relpath).locate()
        assert result == expected
        assert isinstance(result, Path)

    def test_locate_second_prefix(self, two_drives):
        """locate() falls through to second prefix when first has no match."""
        prefix1, prefix2 = two_drives
        relpath = 'sims/run_001/output.log'
        expected = self._make_file(prefix2, relpath)

        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        result = Path(relpath).locate()
        assert result == expected
        assert isinstance(result, Path)

    def test_locate_returns_none(self, two_drives):
        """locate() returns None when no prefix has the file."""
        prefix1, prefix2 = two_drives
        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        assert Path('sims/run_001/missing.log').locate() is None

    def test_locate_absolute_path(self, two_drives):
        """locate() on an absolute path strips the prefix and searches other drives."""
        prefix1, prefix2 = two_drives
        relpath = 'sims/run_001/output.log'
        self._make_file(prefix1, relpath)          # only on drive1
        abs_on_drive2 = Path(prefix2) / relpath    # points at drive2 (doesn't exist)

        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        result = abs_on_drive2.locate()
        assert result is not None
        assert bool(result)

    # ------------------------------------------------------------------
    # locate() — glob patterns and recent
    # ------------------------------------------------------------------

    def test_locate_glob_returns_pattern(self, two_drives):
        """locate() with * returns the pattern (still has *), not an expanded file."""
        prefix1, prefix2 = two_drives
        self._make_file(prefix1, 'sims/run_001/output.log')

        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        result = Path('sims/*/output.log').locate()
        assert result is not None
        assert '*' in str(result)                        # pattern, not expanded
        assert str(prefix1) in str(result)               # correct drive
        assert list(result)                              # expands to real files

    def test_locate_glob_skips_empty_drive(self, two_drives):
        """locate() with * skips drives with no matches and returns the right one."""
        prefix1, prefix2 = two_drives
        self._make_file(prefix2, 'sims/run_001/output.log')  # only on drive2

        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        result = Path('sims/*/output.log').locate()
        assert result is not None
        assert str(prefix2) in str(result)               # fell through to drive2

    def test_locate_recent_newest(self, two_drives):
        """locate(recent=True) returns the most recently modified match."""
        prefix1, prefix2 = two_drives
        self._make_file(prefix1, 'sims/run_001/output.log')
        import time; time.sleep(0.05)
        f2 = self._make_file(prefix2, 'sims/run_001/output.log')

        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        result = Path('sims/run_001/output.log').locate(recent=True)
        assert result == f2

    def test_locate_recent_oldest(self, two_drives):
        """locate(recent=False) returns the oldest (least recently modified) match."""
        prefix1, prefix2 = two_drives
        f1 = self._make_file(prefix1, 'sims/run_001/output.log')
        import time; time.sleep(0.05)
        self._make_file(prefix2, 'sims/run_001/output.log')

        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        result = Path('sims/run_001/output.log').locate(recent=False)
        assert result == f1

    # ------------------------------------------------------------------
    # locate_all() — exact and glob
    # ------------------------------------------------------------------

    def test_locate_all_exact_both_drives(self, two_drives):
        """locate_all() returns one copy per prefix for an exact path."""
        prefix1, prefix2 = two_drives
        relpath = 'sims/run_001/output.log'
        full1 = self._make_file(prefix1, relpath)
        full2 = self._make_file(prefix2, relpath)

        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        results = Path(relpath).locate_all()
        assert len(results) == 2
        assert full1 in results
        assert full2 in results

    def test_locate_all_glob_pattern(self, two_drives):
        """locate_all() collects glob matches across all prefixes."""
        prefix1, prefix2 = two_drives
        f1 = self._make_file(prefix1, 'sims/run_001/output.log')
        f2 = self._make_file(prefix2, 'sims/run_001/output.log')
        f3 = self._make_file(prefix2, 'sims/run_002/output.log')

        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        results = Path('sims/*/output.log').locate_all()
        assert len(results) == 3
        assert f1 in results
        assert f2 in results
        assert f3 in results

    # ------------------------------------------------------------------
    # swap_prefix()
    # ------------------------------------------------------------------

    def test_swap_prefix_by_index(self, two_drives):
        """swap_prefix(int) rewrites the prefix to the indexed entry."""
        prefix1, prefix2 = two_drives
        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        p = Path(prefix2) / 'sims/run_001/output.log'
        result = p.swap_prefix(0)
        assert result == Path(prefix1) / 'sims/run_001/output.log'

    def test_swap_prefix_by_string(self, two_drives):
        """swap_prefix(str) rewrites the prefix to the given string value."""
        prefix1, prefix2 = two_drives
        Path.search_prefixes = [prefix1, prefix2]
        Path.search_common = 'research'

        p = Path(prefix1) / 'sims/run_001/output.log'
        result = p.swap_prefix(str(prefix2))
        assert result == Path(prefix2) / 'sims/run_001/output.log'

    def test_swap_prefix_out_of_range(self, two_drives):
        """swap_prefix() returns None for an out-of-range integer index."""
        prefix1, _ = two_drives
        Path.search_prefixes = [prefix1]
        Path.search_common = 'research'

        p = Path(prefix1) / 'sims/run_001/output.log'
        assert p.swap_prefix(99) is None

    # ------------------------------------------------------------------
    # find_prefix()
    # ------------------------------------------------------------------

    def test_find_prefix_classmethod(self, two_drives):
        """find_prefix() extracts the correct prefix and registers it."""
        prefix1, _ = two_drives
        fake_file = prefix1 / 'sims/run_001/output.log'

        extracted = Path.find_prefix('research', path=fake_file)

        assert extracted == prefix1
        assert Path.search_prefixes == [prefix1]
        assert Path.search_common == 'research'

    def test_find_prefix_accepts_string(self, two_drives):
        """find_prefix() accepts a plain string as the path argument."""
        prefix1, _ = two_drives
        fake_file = str(prefix1 / 'sims/run_001/output.log')

        extracted = Path.find_prefix('research', path=fake_file)

        assert extracted == prefix1

    def test_find_prefix_not_found(self, two_drives):
        """find_prefix() returns None when common dir is not in path."""
        prefix1, _ = two_drives
        result = Path.find_prefix('nonexistent', path=prefix1 / 'sims/output.log')
        assert result is None

    def test_find_prefix_reuses_common(self, two_drives):
        """find_prefix() with common=None reuses the already-set search_common."""
        prefix1, prefix2 = two_drives
        fake1 = prefix1 / 'sims/run_001/output.log'
        fake2 = prefix2 / 'sims/run_001/output.log'

        # First call sets search_common
        Path.find_prefix('research', path=fake1)
        assert Path.search_common == 'research'

        # Second call omits common — should reuse 'research'
        extracted = Path.find_prefix(path=fake2)
        assert extracted == prefix2
        assert len(Path.search_prefixes) == 2

    # ------------------------------------------------------------------
    # Edge cases / configuration errors
    # ------------------------------------------------------------------

    def test_no_prefixes_raises_locate(self):
        """locate() raises RuntimeError when search_prefixes is None."""
        with pytest.raises(RuntimeError, match="search_prefixes"):
            Path('sims/run_001/output.log').locate()

    def test_no_prefixes_raises_locate_all(self):
        """locate_all() raises RuntimeError when search_prefixes is None."""
        with pytest.raises(RuntimeError, match="search_prefixes"):
            Path('sims/run_001/output.log').locate_all()

    def test_empty_prefixes_returns_none(self, two_drives):
        """locate() returns None (not raise) when search_prefixes is [] (configured but empty)."""
        Path.search_prefixes = []
        Path.search_common = 'research'
        assert Path('sims/run_001/output.log').locate() is None
        assert Path('sims/run_001/output.log').locate_all() == []

    def test_absolute_path_no_common_raises(self, two_drives):
        """Operations on an absolute path raise when search_common is None."""
        prefix1, _ = two_drives
        abs_path = Path(prefix1) / 'sims/run_001/output.log'
        Path.search_prefixes = [prefix1]
        # search_common intentionally left as None
        with pytest.raises(RuntimeError, match="search_common"):
            abs_path.locate()

    def test_find_prefix_no_common_raises(self, two_drives):
        """find_prefix() raises RuntimeError when common is None and search_common unset."""
        prefix1, _ = two_drives
        with pytest.raises(RuntimeError, match="search_common"):
            Path.find_prefix(path=prefix1 / 'sims/output.log')
# -*- coding: utf-8 -*-

import pytest
from mooonpy import Path
from mooonpy.molspace.molspace import Molspace
from mooonpy.molspace.atom_styles import Styles
from mooonpy.molspace._files_io.read_lmp_dump import read


class TestReadLmpDump:
    test_dir = Path(__file__).dir() / '../examples/EPON_862/lmp_dump/'
    dump_file = test_dir / 'small_epon_densify.dump'
    data_file = test_dir / 'small_epon.data'

    ## steps in example file
    n_steps = 11
    n_atoms = 468

    @classmethod
    def setup_class(cls):
        """Setup extra fields needed for dump file reading"""
        # Create extra fields for custom dump variables or fixes
        # processor is not supported to test this
        Styles.styles['foo'] = ('proc')
        Styles.defaults.update({'proc': 0})
        Styles.func_aliases.update({'proc': int})

    def test_read_all_steps(self):
        """Test reading all steps from dump file"""
        self.setup_class()
        # Read all steps
        dump_dict = read(self.dump_file, steps=None)

        # Should return a dictionary
        assert isinstance(dump_dict, dict)
        # Should contain multiple timesteps
        assert len(dump_dict) == self.n_steps
        # Each value should be a Molspace
        for step, mol in dump_dict.items():
            assert isinstance(step, int)
            assert hasattr(mol, 'atoms')

    def test_read_single_step(self):
        """Test reading a single specific step"""
        # Read single step (should return Molspace, not dict)
        self.setup_class()
        single_mol = read(self.dump_file, steps=1000)

        # Should return a Molspace directly
        assert hasattr(single_mol, 'atoms')
        assert not isinstance(single_mol, dict)

    def test_read_multiple_steps(self):
        """Test reading multiple specific steps"""
        # Read specific steps
        self.setup_class()
        dump_dict = read(self.dump_file, steps=[0, 1000])

        # Should return dictionary with exactly 2 entries
        assert isinstance(dump_dict, dict)
        assert len(dump_dict) == 2
        assert 0 in dump_dict
        assert 1000 in dump_dict

        # Each should be a Molspace
        for mol in dump_dict.values():
            assert hasattr(mol, 'atoms')

    def test_read_with_base_molspace(self):
        """Test reading dump with a base molspace for inheritance"""
        # Create base molspace with bonded info
        self.setup_class()
        bonded = Molspace(self.data_file)

        # Store original atom type for comparison
        original_atom_type = bonded.atoms[1].type
        assert original_atom_type == 10  # As noted in the example

        # Read dump inheriting from bonded molspace
        dump_dict = read(self.dump_file, steps=[0, 1000], mol=bonded)

        # Should inherit atom types from base molspace
        assert dump_dict[0].atoms[1].type == original_atom_type
        assert dump_dict[1000].atoms[1].type == original_atom_type

    def test_read_with_base_molspace_single_step(self):
        """Test reading single step with base molspace inheritance"""
        # Create base molspace
        self.setup_class()
        bonded = Molspace(self.data_file)
        original_atom_type = bonded.atoms[1].type

        # Read single step with inheritance
        single_mol = read(self.dump_file, steps=1000, mol=bonded)

        # Should inherit atom properties
        assert single_mol.atoms[1].type == original_atom_type
        assert hasattr(single_mol.atoms, 'box')

    def test_box_properties_read(self):
        """Test that box properties are correctly read from dump file"""
        # Read a single step
        self.setup_class()
        mol = read(self.dump_file, steps=0)

        # Check that box properties exist and are non-default
        box = mol.atoms.box

        assert hasattr(box, 'xlo')
        assert hasattr(box, 'xhi')
        assert hasattr(box, 'ylo')
        assert hasattr(box, 'yhi')
        assert hasattr(box, 'zlo')
        assert hasattr(box, 'zhi')

        # Box dimensions should be numeric
        assert isinstance(box.xlo, (int, float))
        assert isinstance(box.xhi, (int, float))
        assert not box.xlo == -0.5
        assert not box.xhi == 0.5
        assert box.xhi > box.xlo  # Sanity check, fails on 0

    def test_atom_coordinates_updated(self):
        """Test that atom coordinates are properly updated from dump file"""
        # Get original coordinates from data file
        self.setup_class()
        bonded = Molspace(self.data_file)
        start_x = bonded.atoms[1].x

        # Read dump file coordinates
        dump_dict = read(self.dump_file, steps=[0, 1000], mol=bonded)

        # Coordinates should be different from original (densification changes positions)
        step0_x = dump_dict[0].atoms[1].x
        step1000_x = dump_dict[1000].atoms[1].x

        # At least one should be different from starting position
        assert step0_x != start_x or step1000_x != start_x

        # Coordinates should exist and be numeric
        assert isinstance(step0_x, (int, float))
        assert isinstance(step1000_x, (int, float))

    def test_nonexistent_step(self):
        """Test reading a step that doesn't exist in the dump file"""
        # Try to read a step that likely doesn't exist
        self.setup_class()
        result = read(self.dump_file, steps=999999)

        # Should return None for non-existent step
        assert result is None

    def test_empty_molspace_base(self):
        """Test reading with None as base molspace (empty)"""

        # Read with None base (should create empty Molspace)
        self.setup_class()
        mol = read(self.dump_file, steps=0, mol=None)

        # Should have atoms from dump but no inherited properties
        assert hasattr(mol, 'atoms')
        assert len(mol.atoms) == self.n_atoms

if __name__ == '__main__':
    pytest.main()
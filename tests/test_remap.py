# -*- coding: utf-8 -*-
"""
Tests for mooonpy.molspace.remap.

Covers:
  * Pure renumbering (offset mode) preserves dict and per-object identities.
  * Round-trip via new2old restores the original state.
  * `contiguous` is a no-op when atoms are already 1..N.
  * Slicing drops un-mapped atoms and any BADI entry referencing them.
  * `reachable` is monotone in depth and self-inclusive.
"""
import pytest
from mooonpy import Molspace, Path


EPON_DIR = Path(__file__).dir() / Path('../examples/EPON_862/atom_typing_Outputs')


def _load_detda():
    """DETDA monomer with atoms+bonds (no FF)."""
    return Molspace(EPON_DIR / Path('detda_typed.data'))


def _load_detda_iff():
    """DETDA with full PCFF-IFF BADI for exercising angle/dihedral/improper remap."""
    return Molspace(EPON_DIR / Path('detda_typed_IFF.data'))


class TestRemapIdentity:
    """Renumbering preserves Python identities of containers and elements."""

    def test_offset_preserves_dict_identity(self):
        mol = _load_detda()
        atoms_dict_id = id(mol.atoms)
        bonds_dict_id = id(mol.bonds)
        first_atom_obj_id = id(mol.atoms[1])

        old2new, _ = mol.build_old2new(mode='offset', offset=100)
        mol.remap_ids(old2new)

        assert id(mol.atoms) == atoms_dict_id
        assert id(mol.bonds) == bonds_dict_id
        assert id(mol.atoms[101]) == first_atom_obj_id
        assert mol.atoms[101].id == 101

    def test_offset_preserves_badi_dict_identity(self):
        mol = _load_detda_iff()
        ids = {
            'angles': id(mol.angles),
            'dihedrals': id(mol.dihedrals),
            'impropers': id(mol.impropers),
        }
        old2new, _ = mol.build_old2new(mode='offset', offset=100)
        mol.remap_ids(old2new)
        assert id(mol.angles) == ids['angles']
        assert id(mol.dihedrals) == ids['dihedrals']
        assert id(mol.impropers) == ids['impropers']


class TestRenumber:
    def test_offset_shifts_all_ids(self):
        mol = _load_detda()
        n_atoms = len(mol.atoms)
        n_bonds = len(mol.bonds)

        old2new, _ = mol.build_old2new(mode='offset', offset=100)
        mol.remap_ids(old2new)

        assert min(mol.atoms) == 101
        assert max(mol.atoms) == 100 + n_atoms
        assert len(mol.bonds) == n_bonds
        for (i, j) in mol.bonds:
            assert i > 100 and j > 100

    def test_round_trip_via_new2old(self):
        mol = _load_detda_iff()
        n_atoms = len(mol.atoms)
        n_bonds = len(mol.bonds)
        n_angles = len(mol.angles)
        n_dihedrals = len(mol.dihedrals)
        n_impropers = len(mol.impropers)

        old2new, new2old = mol.build_old2new(mode='offset', offset=500)
        mol.remap_ids(old2new)
        mol.remap_ids(new2old)

        assert len(mol.atoms) == n_atoms
        assert len(mol.bonds) == n_bonds
        assert len(mol.angles) == n_angles
        assert len(mol.dihedrals) == n_dihedrals
        assert len(mol.impropers) == n_impropers
        assert min(mol.atoms) == 1
        assert max(mol.atoms) == n_atoms

    def test_contiguous_is_noop_when_already_contiguous(self):
        mol = _load_detda()
        n = len(mol.atoms)
        mol.contiguous()
        assert set(mol.atoms.keys()) == set(range(1, n + 1))

    def test_contiguous_after_offset_returns_to_one_based(self):
        mol = _load_detda()
        n = len(mol.atoms)
        old2new, _ = mol.build_old2new(mode='offset', offset=1000)
        mol.remap_ids(old2new)
        mol.contiguous()
        assert set(mol.atoms.keys()) == set(range(1, n + 1))

    def test_badi_ordered_field_remapped(self):
        mol = _load_detda_iff()
        old2new, _ = mol.build_old2new(mode='offset', offset=100)
        mol.remap_ids(old2new)
        for angle in mol.angles.values():
            assert all(v > 100 for v in angle.ordered)
        for dihedral in mol.dihedrals.values():
            assert all(v > 100 for v in dihedral.ordered)
        for improper in mol.impropers.values():
            assert all(v > 100 for v in improper.ordered)


class TestSlice:
    def test_slice_identity_drops_atoms_and_bonds(self):
        mol = _load_detda()
        keep = set(range(1, 16))
        mol.slice(keep, renumber=False)

        assert set(mol.atoms.keys()) == keep
        for (i, j) in mol.bonds:
            assert i in keep and j in keep

    def test_slice_contiguous_renumbers_to_one_based(self):
        mol = _load_detda()
        keep = set(range(1, 16))
        mol.slice(keep, renumber=True)
        assert set(mol.atoms.keys()) == set(range(1, 16))

    def test_slice_cascades_through_all_badi(self):
        mol = _load_detda_iff()
        keep = set(range(1, 11))
        mol.slice(keep, renumber=False)
        for key in mol.bonds:
            assert all(v in keep for v in key)
        for key in mol.angles:
            assert all(v in keep for v in key)
        for key in mol.dihedrals:
            assert all(v in keep for v in key)
        for key in mol.impropers:
            assert all(v in keep for v in key)

    def test_slice_preserves_atom_object_identity(self):
        mol = _load_detda()
        original = {id_: id(mol.atoms[id_]) for id_ in range(1, 11)}
        mol.slice(set(range(1, 11)), renumber=False)
        for id_, obj_id in original.items():
            assert id(mol.atoms[id_]) == obj_id


class TestReachable:
    def test_self_inclusive(self):
        mol = _load_detda()
        assert 1 in mol.reachable(1, depth=1)
        assert 1 in mol.reachable(1, depth=0)

    def test_depth_zero_is_just_start(self):
        mol = _load_detda()
        assert mol.reachable(1, depth=0) == {1}

    def test_monotone_in_depth(self):
        mol = _load_detda()
        r1 = mol.reachable(1, depth=1)
        r2 = mol.reachable(1, depth=2)
        r3 = mol.reachable(1, depth=3)
        assert r1.issubset(r2)
        assert r2.issubset(r3)

    def test_negative_depth_raises(self):
        mol = _load_detda()
        with pytest.raises(ValueError):
            mol.reachable(1, depth=-1)

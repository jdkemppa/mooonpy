# -*- coding: utf-8 -*-
"""
Tests for Molspace.combine.

Combine appends another Molspace into ``self`` with id shifting (no collisions)
and an optional translation.
"""
import pytest
from mooonpy import Molspace, Path


EPON_DIR = Path(__file__).dir() / Path('../examples/EPON_862/atom_typing_Outputs')


def _load_detda_iff():
    return Molspace(EPON_DIR / Path('detda_typed_IFF.data'))


class TestCombine:
    def test_doubles_atoms_and_bonds(self):
        a = _load_detda_iff()
        b = _load_detda_iff()
        n_atoms_a = len(a.atoms)
        n_bonds_a = len(a.bonds)
        n_angles_a = len(a.angles)
        n_dihedrals_a = len(a.dihedrals)

        a.combine(b, vect=(10.0, 0.0, 0.0))

        assert len(a.atoms) == 2 * n_atoms_a
        assert len(a.bonds) == 2 * n_bonds_a
        assert len(a.angles) == 2 * n_angles_a
        assert len(a.dihedrals) == 2 * n_dihedrals_a

    def test_no_id_collisions(self):
        a = _load_detda_iff()
        b = _load_detda_iff()
        n = len(a.atoms)
        a.combine(b, vect=(10.0, 0.0, 0.0))
        assert set(a.atoms.keys()) == set(range(1, 2 * n + 1))

    def test_other_not_mutated(self):
        a = _load_detda_iff()
        b = _load_detda_iff()
        n_b_before = len(b.atoms)
        first_b_x = b.atoms[1].x
        a.combine(b, vect=(10.0, 0.0, 0.0))
        assert len(b.atoms) == n_b_before
        assert b.atoms[1].x == first_b_x  # b's coordinates untouched

    def test_translation_applied_to_appended_copy(self):
        a = _load_detda_iff()
        b = _load_detda_iff()
        n = len(a.atoms)
        original_xs = [b.atoms[id_].x for id_ in b.atoms]

        a.combine(b, vect=(10.0, 0.0, 0.0), move_mode='offset')

        # Appended atoms have ids n+1 .. 2n; their x should be original + 10.
        for k, original_x in enumerate(original_xs, start=1):
            assert a.atoms[n + k].x == pytest.approx(original_x + 10.0)

    def test_no_vect_keeps_appended_coords(self):
        a = _load_detda_iff()
        b = _load_detda_iff()
        n = len(a.atoms)
        original_xs = [b.atoms[id_].x for id_ in b.atoms]

        a.combine(b, vect=None)

        for k, original_x in enumerate(original_xs, start=1):
            assert a.atoms[n + k].x == pytest.approx(original_x)

    def test_two_clusters_after_combine(self):
        a = _load_detda_iff()
        b = _load_detda_iff()
        a.combine(b, vect=(50.0, 0.0, 0.0))
        # Two disconnected monomers → at least two clusters.
        assert len(a.clusters) == 2

    def test_combine_into_empty(self):
        a = Molspace()
        b = _load_detda_iff()
        n_b = len(b.atoms)
        shift = a.combine(b, vect=None)
        assert shift == 0
        assert len(a.atoms) == n_b
        assert set(a.atoms.keys()) == set(range(1, n_b + 1))

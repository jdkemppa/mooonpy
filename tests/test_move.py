# -*- coding: utf-8 -*-
"""
Tests for Atoms.move / Atoms.centroid / Atoms.COM and box-aware scale.
"""
import pytest
from mooonpy import Molspace, Path


EPON_DIR = Path(__file__).dir() / Path('../examples/EPON_862/atom_typing_Outputs')


def _load_detda():
    return Molspace(EPON_DIR / Path('detda_typed.data'))


def _load_detda_iff():
    """Has masses; usable for COM."""
    return Molspace(EPON_DIR / Path('detda_typed_IFF.data'))


class TestMove:
    def test_none_is_noop(self):
        mol = _load_detda()
        before = [(a.x, a.y, a.z) for a in mol.atoms.values()]
        mol.atoms.move(None)
        after = [(a.x, a.y, a.z) for a in mol.atoms.values()]
        assert before == after

    def test_offset_shifts_all_atoms(self):
        mol = _load_detda()
        before = [(a.x, a.y, a.z) for a in mol.atoms.values()]
        mol.atoms.move((1.0, 2.0, 3.0), mode='offset')
        after = [(a.x, a.y, a.z) for a in mol.atoms.values()]
        for b, a in zip(before, after):
            assert a[0] == pytest.approx(b[0] + 1.0)
            assert a[1] == pytest.approx(b[1] + 2.0)
            assert a[2] == pytest.approx(b[2] + 3.0)

    def test_offset_does_not_scale_box(self):
        mol = _load_detda()
        box_before = (mol.atoms.box.xlo, mol.atoms.box.xhi)
        mol.atoms.move((10.0, 0.0, 0.0), mode='offset')
        assert (mol.atoms.box.xlo, mol.atoms.box.xhi) == box_before

    def test_centroid_lands_at_target(self):
        mol = _load_detda()
        mol.atoms.move((5.0, -3.0, 2.5), mode='centroid')
        cx, cy, cz = mol.atoms.centroid()
        assert cx == pytest.approx(5.0, abs=1e-10)
        assert cy == pytest.approx(-3.0, abs=1e-10)
        assert cz == pytest.approx(2.5, abs=1e-10)

    def test_scale_multiplies_each_axis_and_box(self):
        mol = _load_detda()
        before = [(a.x, a.y, a.z) for a in mol.atoms.values()]
        box_before = (mol.atoms.box.xlo, mol.atoms.box.xhi,
                      mol.atoms.box.ylo, mol.atoms.box.yhi,
                      mol.atoms.box.zlo, mol.atoms.box.zhi)
        mol.atoms.move((2.0, 1.0, 0.5), mode='scale')
        after = [(a.x, a.y, a.z) for a in mol.atoms.values()]
        for b, a in zip(before, after):
            assert a[0] == pytest.approx(b[0] * 2.0)
            assert a[1] == pytest.approx(b[1] * 1.0)
            assert a[2] == pytest.approx(b[2] * 0.5)
        # Box scales about origin.
        assert mol.atoms.box.xlo == pytest.approx(box_before[0] * 2.0)
        assert mol.atoms.box.xhi == pytest.approx(box_before[1] * 2.0)
        assert mol.atoms.box.ylo == pytest.approx(box_before[2] * 1.0)
        assert mol.atoms.box.yhi == pytest.approx(box_before[3] * 1.0)
        assert mol.atoms.box.zlo == pytest.approx(box_before[4] * 0.5)
        assert mol.atoms.box.zhi == pytest.approx(box_before[5] * 0.5)

    def test_atom_ids_subset(self):
        """``atom_ids`` restricts the move to a subset of atoms."""
        mol = _load_detda()
        before = {id_: (a.x, a.y, a.z) for id_, a in mol.atoms.items()}
        mol.atoms.move((1.0, 0.0, 0.0), mode='offset', atom_ids=[1, 2])
        for id_, a in mol.atoms.items():
            if id_ in (1, 2):
                assert a.x == pytest.approx(before[id_][0] + 1.0)
            else:
                assert a.x == pytest.approx(before[id_][0])

    def test_unknown_mode_raises(self):
        mol = _load_detda()
        with pytest.raises(ValueError):
            mol.atoms.move((0, 0, 0), mode='wat')


class TestCentroid:
    def test_centroid_matches_manual_computation(self):
        mol = _load_detda()
        n = len(mol.atoms)
        cx_manual = sum(a.x for a in mol.atoms.values()) / n
        cy_manual = sum(a.y for a in mol.atoms.values()) / n
        cz_manual = sum(a.z for a in mol.atoms.values()) / n
        cx, cy, cz = mol.atoms.centroid()
        assert cx == pytest.approx(cx_manual)
        assert cy == pytest.approx(cy_manual)
        assert cz == pytest.approx(cz_manual)

    def test_centroid_subset(self):
        mol = _load_detda()
        c1, _, _ = mol.atoms.centroid(atom_ids=[1])
        assert c1 == pytest.approx(mol.atoms[1].x)


class TestCOM:
    def test_com_with_ff_masses(self):
        mol = _load_detda_iff()
        assert len(mol.ff.masses) > 0
        cx, cy, cz = mol.atoms.COM(mol.ff.masses)
        # Mass-weighted manual computation
        masses = {t: p.coeffs[0] for t, p in mol.ff.masses.items()}
        total = sum(masses[a.type] for a in mol.atoms.values())
        mx = sum(masses[a.type] * a.x for a in mol.atoms.values()) / total
        my = sum(masses[a.type] * a.y for a in mol.atoms.values()) / total
        mz = sum(masses[a.type] * a.z for a in mol.atoms.values()) / total
        assert cx == pytest.approx(mx)
        assert cy == pytest.approx(my)
        assert cz == pytest.approx(mz)

    def test_com_with_plain_dict(self):
        mol = _load_detda_iff()
        masses = {t: p.coeffs[0] for t, p in mol.ff.masses.items()}
        ff_com = mol.atoms.COM(mol.ff.masses)
        plain_com = mol.atoms.COM(masses)
        for a, b in zip(ff_com, plain_com):
            assert a == pytest.approx(b)

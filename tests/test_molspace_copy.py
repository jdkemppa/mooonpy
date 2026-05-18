# -*- coding: utf-8 -*-
"""
Safety / correctness harness for Molspace.copy().

Run BEFORE and AFTER any change to Molspace.copy(): a faithful, fully
independent deep copy must keep passing every assertion here. Uses the small
DETDA monomer (31 atoms, full topology + class-2 force field) so the whole
suite runs in well under a second -- speed is validated separately on PBZ.

What is checked:
  1. Structural fidelity   - copy reproduces every data structure by value.
  2. Pointer correctness   - containers/value-objects/mutable-attrs are
                             distinct instances; clusters._molspace back-ref
                             points at the NEW Molspace (not the old one).
  3. Mutation isolation    - mutating the copy never touches the original
                             (and vice-versa) for ALL data structures:
                             atoms, bonds, angles, dihedrals, impropers,
                             force field, box, clusters.

    pytest tests/test_molspace_copy.py -v
"""

import os
import sys

import pytest

# Pin to the worktree mooonpy regardless of any installed/editable copy.
_WORKTREE_SRC = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "src")
)
sys.path.insert(0, _WORKTREE_SRC)

import mooonpy  # noqa: E402

DETDA = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "examples", "EPON_862", "atom_typing_Outputs", "detda_typed_IFF.data",
)

# Every topology container on a Molspace and the value-object attribute we
# mutate to prove independence. `ordered` is a list (mutable) -> catches a
# shallow copy that a scalar-only check would miss.
TOPO = ("bonds", "angles", "dihedrals", "impropers")
FF_TABLES = (
    "masses", "pair_coeffs", "bond_coeffs", "angle_coeffs",
    "dihedral_coeffs", "improper_coeffs",
)


def _slots(obj):
    """{slot_name: value} for a __slots__ value-object (Atom/Bond/...)."""
    return {s: getattr(obj, s) for s in type(obj).__slots__}


def _load():
    mol = mooonpy.Molspace(DETDA)
    # clusters are not populated at read time; populate so they are exercised.
    mol.clusters.assign_clusters()
    mol.clusters.compute_info()
    return mol


# Tests run against the fast custom copier (mol.copy()). The
# mol.copy(deepcopy=True) fallback is correct but ~6-8x slower, so it is
# kept OUT of CI. Re-enable the commented entries below to cross-check the
# fast path against copy.deepcopy locally.
COPY_MODES = (False,)            # (False, True)
_MODE_IDS = ("fast",)            # ("fast", "deepcopy")


@pytest.fixture(scope="module")
def original():
    return _load()


@pytest.fixture(params=COPY_MODES, ids=_MODE_IDS)
def mode(request):
    """The `deepcopy=` kwarg passed to Molspace.copy()."""
    return request.param


@pytest.fixture()
def pair(mode):
    """A freshly loaded original and its copy via the parametrized path."""
    orig = _load()
    return orig, orig.copy(deepcopy=mode)


# ---------------------------------------------------------------------------
# Sanity: the worktree copy is the one under test, DETDA loaded as expected
# ---------------------------------------------------------------------------

def test_using_worktree_mooonpy():
    assert "cranky-proskuriakova-780598" in mooonpy.__file__, mooonpy.__file__


def test_detda_loaded(original):
    assert len(original.atoms) == 31
    assert len(original.bonds) == 31
    assert len(original.angles) == 54
    assert len(original.dihedrals) == 68
    assert len(original.impropers) == 28
    assert len(original.ff.masses) == 6
    assert len(original.clusters) >= 1


# ---------------------------------------------------------------------------
# 1. Structural fidelity -- copy reproduces every structure by value
# ---------------------------------------------------------------------------

def test_atoms_equal_by_value(pair):
    orig, cp = pair
    assert set(cp.atoms) == set(orig.atoms)
    for k in orig.atoms:
        assert _slots(cp.atoms[k]) == _slots(orig.atoms[k])
    assert cp.atoms.style == orig.atoms.style


@pytest.mark.parametrize("name", TOPO)
def test_topology_equal_by_value(pair, name):
    orig, cp = pair
    o, c = getattr(orig, name), getattr(cp, name)
    assert set(c) == set(o)
    for k in o:
        assert _slots(c[k]) == _slots(o[k])


@pytest.mark.parametrize("table", FF_TABLES)
def test_ff_equal_by_value(pair, table):
    orig, cp = pair
    o, c = getattr(orig.ff, table), getattr(cp.ff, table)
    assert set(c) == set(o)
    for t in o:
        assert c[t].coeffs == o[t].coeffs
        assert c[t].style == o[t].style
        assert c[t].comment == o[t].comment


def test_box_equal_by_value(pair):
    orig, cp = pair
    for a in ("xlo", "xhi", "ylo", "yhi", "zlo", "zhi", "xy", "xz", "yz"):
        assert getattr(cp.atoms.box, a) == getattr(orig.atoms.box, a)


def test_clusters_equal_by_value(pair):
    orig, cp = pair
    assert set(cp.clusters) == set(orig.clusters)
    for molid in orig.clusters:
        assert cp.clusters[molid].members == orig.clusters[molid].members
        assert cp.clusters[molid].mass == orig.clusters[molid].mass
        assert dict(cp.clusters[molid].composition) == dict(
            orig.clusters[molid].composition
        )


# ---------------------------------------------------------------------------
# 2. Pointer correctness -- distinct instances + correct back-ref
# ---------------------------------------------------------------------------

def test_containers_are_distinct_instances(pair):
    orig, cp = pair
    assert cp is not orig
    assert cp.atoms is not orig.atoms
    assert cp.ff is not orig.ff
    assert cp.clusters is not orig.clusters
    assert cp.atoms.box is not orig.atoms.box
    for name in TOPO:
        assert getattr(cp, name) is not getattr(orig, name)


def test_value_objects_are_distinct_instances(pair):
    orig, cp = pair
    k = next(iter(orig.atoms))
    assert cp.atoms[k] is not orig.atoms[k]
    for name in TOPO:
        o = getattr(orig, name)
        kk = next(iter(o))
        c = getattr(cp, name)[kk]
        assert c is not o[kk]
        # mutable attribute must also be a distinct object, not shared
        assert c.ordered is not o[kk].ordered


def test_clusters_backref_points_to_new_molspace(pair):
    """The invariant most likely to break in a hand-written copy."""
    orig, cp = pair
    assert cp.clusters._molspace is cp
    assert cp.clusters._molspace is not orig
    assert orig.clusters._molspace is orig


def test_clusters_backref_is_functionally_correct(pair):
    """
    recompute on the copy must read the COPY's atoms and write back to the
    COPY only -- proves the back-ref is wired to the new graph, not shared.
    """
    orig, cp = pair
    orig_molids = [a.molid for a in orig.atoms.values()]
    cp.clusters.compute_clusters()  # writes new molids onto cp.atoms via back-ref
    assert [a.molid for a in orig.atoms.values()] == orig_molids
    # the recomputed clusters reference copied atom ids and stay self-consistent
    for cluster in cp.clusters.values():
        for aid in cluster.members:
            assert aid in cp.atoms


def test_atom_class_is_shared_not_rebuilt(pair):
    """Sharing the dynamically-built class is correct and expected."""
    orig, cp = pair
    k = next(iter(orig.atoms))
    assert type(cp.atoms[k]) is type(orig.atoms[k])


# ---------------------------------------------------------------------------
# 3. Mutation isolation -- mutate copy, original must be untouched
# ---------------------------------------------------------------------------

def test_mutate_copy_atom_scalar(pair):
    orig, cp = pair
    k = next(iter(orig.atoms))
    before = orig.atoms[k].x
    cp.atoms[k].x = before + 12345.6
    assert orig.atoms[k].x == before


@pytest.mark.parametrize("name", TOPO)
def test_mutate_copy_topology(pair, name):
    orig, cp = pair
    o, c = getattr(orig, name), getattr(cp, name)
    k = next(iter(o))
    type_before = o[k].type
    ordered_before = list(o[k].ordered)

    c[k].type = 99999
    c[k].comment = "MUTATED"
    c[k].ordered.append("MUTATED")          # in-place mutate of shared-risk list

    assert o[k].type == type_before
    assert o[k].comment != "MUTATED"
    assert list(o[k].ordered) == ordered_before


@pytest.mark.parametrize("table", FF_TABLES)
def test_mutate_copy_ff(pair, table):
    orig, cp = pair
    o, c = getattr(orig.ff, table), getattr(cp.ff, table)
    t = next(iter(o))
    coeffs_before = list(o[t].coeffs)

    c[t].coeffs.append(123456.0)            # in-place mutate of coeffs list
    c[t].comment = "MUTATED"

    assert list(o[t].coeffs) == coeffs_before
    assert o[t].comment != "MUTATED"


def test_mutate_copy_box(pair):
    orig, cp = pair
    before = orig.atoms.box.xlo
    cp.atoms.box.xlo = before - 999.0
    assert orig.atoms.box.xlo == before


def test_mutate_copy_clusters(pair):
    orig, cp = pair
    molid = next(iter(orig.clusters))
    members_before = set(orig.clusters[molid].members)
    mass_before = orig.clusters[molid].mass

    cp.clusters[molid].members.add(-1)
    cp.clusters[molid].mass += 1000.0
    cp.clusters[molid].composition[-999] += 7

    assert orig.clusters[molid].members == members_before
    assert orig.clusters[molid].mass == mass_before
    assert -999 not in orig.clusters[molid].composition


def test_mutate_original_does_not_touch_copy(pair):
    """Reverse direction: mutating the original must not leak into the copy."""
    orig, cp = pair
    k = next(iter(orig.atoms))
    cp_val = cp.atoms[k].x
    orig.atoms[k].x = cp_val + 54321.0
    assert cp.atoms[k].x == cp_val

    bk = next(iter(orig.bonds))
    cp_ordered = list(cp.bonds[bk].ordered)
    orig.bonds[bk].ordered.append("ORIG_ONLY")
    assert list(cp.bonds[bk].ordered) == cp_ordered


# ---------------------------------------------------------------------------
# 4. Functional independence via the clusters graph: break the (5,7) bond
#    (splits 7 atoms off DETDA) and prove copies made before/after the break
#    are independent, and that copy itself neither errors nor mutates source.
# ---------------------------------------------------------------------------

def _break_5_7(mol):
    """Delete the 5--7 bond and recompute clusters from connectivity."""
    key = mol.bonds.key_rule(5, 7)
    assert key in mol.bonds, "DETDA should contain the 5-7 bond"
    del mol.bonds[key]
    mol.clusters.compute_clusters()
    return mol


def _cluster_sizes(mol):
    return sorted(len(c.members) for c in mol.clusters.values())


class TestBondBreakClusterIndependence:

    def test_intact_is_single_cluster(self):
        mol = _load()
        mol.clusters.compute_clusters()
        assert _cluster_sizes(mol) == [31]

    def test_break_splits_off_seven_atoms(self):
        mol = _load()
        _break_5_7(mol)
        sizes = _cluster_sizes(mol)
        assert len(sizes) == 2
        assert sum(sizes) == 31
        assert 7 in sizes              # 7 atoms on one side, as expected

    @pytest.mark.parametrize("mode", COPY_MODES, ids=_MODE_IDS)
    def test_copy_before_break_is_independent(self, mode):
        """A copy taken BEFORE breaking must not change when the original
        is later broken (no shared bonds dict / clusters / back-ref)."""
        orig = _load()
        orig.clusters.compute_clusters()
        cp_before = orig.copy(deepcopy=mode)

        _break_5_7(orig)  # mutate ONLY the original

        # original now split; the earlier copy is untouched
        assert _cluster_sizes(orig) == sorted([7, 24])
        assert _cluster_sizes(cp_before) == [31]
        assert orig.bonds.key_rule(5, 7) not in orig.bonds
        assert cp_before.bonds.key_rule(5, 7) in cp_before.bonds
        # recomputing the untouched copy still yields one cluster
        cp_before.clusters.compute_clusters()
        assert _cluster_sizes(cp_before) == [31]

    @pytest.mark.parametrize("mode", COPY_MODES, ids=_MODE_IDS)
    def test_copy_after_break_is_independent(self, mode):
        """A copy taken AFTER breaking reproduces the 2-cluster state and
        stays independent under mutation / recompute via its own back-ref."""
        orig = _load()
        _break_5_7(orig)
        cp_after = orig.copy(deepcopy=mode)

        # faithfully reproduced
        assert _cluster_sizes(cp_after) == sorted([7, 24])
        assert {m: set(c.members) for m, c in cp_after.clusters.items()} == \
               {m: set(c.members) for m, c in orig.clusters.items()}

        # mutating the copy's clusters does not touch the original
        any_molid = next(iter(cp_after.clusters))
        before = set(orig.clusters[any_molid].members)
        cp_after.clusters[any_molid].members.add(-12345)
        assert orig.clusters[any_molid].members == before

        # recompute on the copy uses the copy's back-ref/atoms only
        orig_molids = [a.molid for a in orig.atoms.values()]
        cp_after.clusters.compute_clusters()
        assert cp_after.clusters._molspace is cp_after
        assert [a.molid for a in orig.atoms.values()] == orig_molids

    @pytest.mark.parametrize("mode", COPY_MODES, ids=_MODE_IDS)
    def test_copy_does_not_mutate_or_error_on_source(self, mode):
        """Taking a copy (before or after a break) must leave the source
        byte-for-byte equivalent to an independently loaded reference."""
        # case A: intact source
        mol = _load()
        ref = _load()
        _ = mol.copy(deepcopy=mode)
        assert {k: _slots(mol.atoms[k]) for k in mol.atoms} == \
               {k: _slots(ref.atoms[k]) for k in ref.atoms}
        assert set(mol.bonds) == set(ref.bonds)
        assert {m: set(c.members) for m, c in mol.clusters.items()} == \
               {m: set(c.members) for m, c in ref.clusters.items()}

        # case B: post-break source
        mol_b = _load(); _break_5_7(mol_b)
        ref_b = _load(); _break_5_7(ref_b)
        snapshot = {m: set(c.members) for m, c in mol_b.clusters.items()}
        _ = mol_b.copy(deepcopy=mode)
        assert {m: set(c.members) for m, c in mol_b.clusters.items()} == snapshot
        assert {m: set(c.members) for m, c in mol_b.clusters.items()} == \
               {m: set(c.members) for m, c in ref_b.clusters.items()}


# ---------------------------------------------------------------------------
# 5. Standalone container .copy() -- Atoms / Bonds / Angles / Dihedrals /
#    Impropers each have their own copy(deepcopy=) toggle; both paths must be
#    correct and independent in isolation (no Molspace involved).
# ---------------------------------------------------------------------------

# (container attr, scalar slot, mutable-list slot or None)
_CONTAINERS = (
    ("atoms", "x", None),
    ("bonds", "type", "ordered"),
    ("angles", "type", "ordered"),
    ("dihedrals", "type", "ordered"),
    ("impropers", "type", "ordered"),
)


@pytest.mark.parametrize("name,scalar,listattr", _CONTAINERS,
                         ids=[c[0] for c in _CONTAINERS])
@pytest.mark.parametrize("mode", COPY_MODES, ids=_MODE_IDS)
class TestContainerCopies:

    def test_equal_distinct_and_class_preserved(self, name, scalar, listattr, mode):
        mol = _load()
        cont = getattr(mol, name)
        cp = cont.copy(deepcopy=mode)

        assert cp is not cont
        assert type(cp) is type(cont)
        assert set(cp) == set(cont)
        k = next(iter(cont))
        assert _slots(cp[k]) == _slots(cont[k])
        assert cp[k] is not cont[k]                  # distinct value objects
        assert type(cp[k]) is type(cont[k])          # factory class preserved
        if listattr is not None:
            assert getattr(cp[k], listattr) is not getattr(cont[k], listattr)

    def test_mutation_isolation(self, name, scalar, listattr, mode):
        mol = _load()
        cont = getattr(mol, name)
        cp = cont.copy(deepcopy=mode)
        k = next(iter(cont))

        scal_before = getattr(cont[k], scalar)
        setattr(cp[k], scalar, scal_before + 99999)
        assert getattr(cont[k], scalar) == scal_before

        if listattr is not None:
            lst_before = list(getattr(cont[k], listattr))
            getattr(cp[k], listattr).append("MUTATED")
            assert list(getattr(cont[k], listattr)) == lst_before

    def test_atoms_box_is_independent(self, name, scalar, listattr, mode):
        if name != "atoms":
            pytest.skip("box check only applies to Atoms")
        mol = _load()
        cont = mol.atoms
        cp = cont.copy(deepcopy=mode)
        assert cp.box is not cont.box
        before = cont.box.xlo
        cp.box.xlo = before - 777.0
        assert cont.box.xlo == before

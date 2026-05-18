# -*- coding: utf-8 -*-
"""
Tests for Fragment + ReactionTemplate against the actual EPON 862 chemistry.

The two crosslinking reactions:

    Rxn 1 — primary amine + epoxide -> secondary amine + alcohol

        N(H)(H)R + C(H)(R')-O-C(H)(R'')  ->  N(H)(R)-C(H)(R')-C(OH)(H)(R'')

      Mechanistically, the amine N attacks one of the epoxide ring carbons
      and the epoxide ring opens at the attacked C-O bond. One amine H
      transfers from N to the former epoxide O, completing the alcohol.
      All atom ids are preserved between pre and post; the H that transfers
      keeps its id but appears bonded to a different heavy atom.

    Rxn 2 — secondary amine + epoxide -> tertiary amine + alcohol

        Same mechanism on the now-secondary amine produced by Rxn 1. We
        construct the input for Rxn 2 by running Rxn 1 first and feeding its
        ``post`` Fragment straight into a second template build.

``ReactionTemplate.pre`` and ``.post`` are :class:`Fragment` objects (access
the Molspace via ``.mol``); they carry the combined edge set and the severed-
neighbor NTAs so a follow-up reaction can slice the post directly.

Modern atom_typing-output NTAs (PCFF/PCFF-IFF naming):
    - DETDA primary amine N: 'nn'  (kept through every reaction here; this FF
      practice does not tag secondary/tertiary amines separately)
    - DETDA amine H:         'hn2'
    - DGEBF epoxide O:       'o3e'
    - DGEBF epoxide C-H:     'c3h'
    - DGEBF methylene C:     'c2'
    - aliphatic H on C:      'hpan'
    - alcohol O:             'oh' (post-reaction)
    - alcohol/hydroxyl H:    'ho' (post-reaction)
    - opened epoxide C:      'c2' (attacked) / 'c1' (other)
"""
import os
import tempfile
import warnings
import pytest
from mooonpy import Molspace, Path, ReactionTemplate
from mooonpy.template import Fragment


EPON_DIR = Path(__file__).dir() / Path('../examples/EPON_862/atom_typing_Outputs')


def _load_detda():
    return Molspace(EPON_DIR / Path('detda_typed.data'))


def _load_dgebf():
    return Molspace(EPON_DIR / Path('dgebf_typed.data'))


def _nta(atom):
    """Pull the NTA token off the verbose comment field that atom_typing writes."""
    if not atom.comment:
        return ''
    parts = atom.comment.split()
    if 'type:' in parts:
        return parts[parts.index('type:') + 1]
    return parts[-1]


def _amine_hs(mol, nitrogen_id):
    """User gives the N id; walk one bond to find the amine H atoms."""
    return [nbr for nbr in mol.generate_graph()[nitrogen_id]
            if _nta(mol.atoms[nbr]) == 'hn2']


def _epoxide_partners(mol, attacked_c_id):
    """User gives the attacked epoxide C; walk one bond to find the O and the
    other epoxide C."""
    nbrs = mol.generate_graph()[attacked_c_id]
    o_id = next(nbr for nbr in nbrs if _nta(mol.atoms[nbr]) == 'o3e')
    other_c = next(nbr for nbr in nbrs
                   if _nta(mol.atoms[nbr]) == 'c3h' and nbr != attacked_c_id)
    return o_id, other_c


# Modern atom_typing IDs in detda_typed.data / dgebf_typed.data.
DETDA_N_PRIMARY = 12     # one of two primary amines; the other is N=13
# The terminal CH2 epoxide carbon (two H neighbors, ids 24 & 25). Its ring
# partner C=3 is the internal CH (one H, plus the methylene C=4).
DGEBF_EPOXIDE_C = 2      # the C the amine attacks

# Edges sit ≥3 bonds out from the reactive atom so every dihedral that crosses
# the new/broken bond has all four atoms inside the template (4 atoms = 3
# bonds), which fix bond/react needs to match parameters. DGEBF gets one extra
# step to cleanly reach the phenyl ring across the methylene-O linker.
DETDA_DEPTH = 3
DGEBF_DEPTH = 4


class TestGraphHelpers:
    """User gives the N or epoxide-C id; helpers find the H/O via graph."""

    def test_amine_hs_discovered_from_n(self):
        detda = _load_detda()
        hs = sorted(_amine_hs(detda, DETDA_N_PRIMARY))
        assert hs == [28, 29]

    def test_epoxide_partners_discovered_from_c(self):
        dgebf = _load_dgebf()
        o_id, other_c = _epoxide_partners(dgebf, DGEBF_EPOXIDE_C)
        assert o_id == 1
        assert other_c == 3


class TestFragmentChemistry:
    """The fragment must enclose the reactive atoms plus enough atoms beyond
    the reactive bond for dihedral matching, with single-interior-bond edges."""

    def test_detda_amine_fragment_depth_3(self):
        detda = _load_detda()
        frag = Fragment.from_initiator_steps(detda, DETDA_N_PRIMARY, n_steps=DETDA_DEPTH)
        # Reactive atoms inside.
        assert DETDA_N_PRIMARY in frag.mol.atoms
        for h in (28, 29):
            assert h in frag.mol.atoms
        # Every edge has exactly one bond into the interior.
        graph = frag.mol.generate_graph()
        for edge in frag.edges:
            n_interior = sum(1 for n in graph[edge] if n in frag.mol.atoms)
            assert n_interior == 1, f"edge {edge} has {n_interior} interior bonds"

    def test_dgebf_epoxide_fragment_depth_4(self):
        dgebf = _load_dgebf()
        frag = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)
        # Epoxide ring (O=1, C=2, C=3) and the H on each epoxide C are inside.
        for atom_id in (1, 2, 3, 24, 25, 26):
            assert atom_id in frag.mol.atoms
        # Phenyl-side reach: the ether O=5 and ring C=6 must both be inside,
        # so dihedrals across the C(epox)-C(methylene)-O-C(ring) linker fit.
        for atom_id in (5, 6):
            assert atom_id in frag.mol.atoms
        graph = frag.mol.generate_graph()
        for edge in frag.edges:
            n_interior = sum(1 for n in graph[edge] if n in frag.mol.atoms)
            assert n_interior == 1, f"edge {edge} has {n_interior} interior bonds"


class TestReaction1AmineEpoxide:
    """Primary amine + epoxide -> secondary amine + alcohol, with H-transfer."""

    def _build(self):
        detda = _load_detda()
        dgebf = _load_dgebf()

        amine_hs = _amine_hs(detda, DETDA_N_PRIMARY)
        h_transfer = amine_hs[0]  # arbitrarily pick one H to transfer
        # The unchanged amine H stays bonded to N — useful for assertions.
        h_remaining = amine_hs[1]
        epox_o, other_epox_c = _epoxide_partners(dgebf, DGEBF_EPOXIDE_C)

        A = Fragment.from_initiator_steps(detda, DETDA_N_PRIMARY, n_steps=DETDA_DEPTH)
        B = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)

        rxn = ReactionTemplate(
            A, B,
            A_name='detda', B_name='dgebf',
            breaks=[
                ('detda', DETDA_N_PRIMARY, h_transfer),  # break the N-H of the H that transfers
                ('dgebf', epox_o, DGEBF_EPOXIDE_C),      # break the epoxide O-C on the attacked side
            ],
            makes=[
                ('detda', DETDA_N_PRIMARY, 'dgebf', DGEBF_EPOXIDE_C),  # new N-C
                ('detda', h_transfer, 'dgebf', epox_o),                # H transfers from N to O
            ],
            changed_typ=[
                # N keeps NTA 'nn' through both reactions in this FF practice
                # (we don't tag primary/secondary/tertiary amines separately).
                ('detda', h_transfer, 'ho'),          # amine H -> hydroxyl H
                ('dgebf', epox_o, 'oh'),              # epoxide O -> hydroxyl O
                ('dgebf', DGEBF_EPOXIDE_C, 'c2'),     # opened epoxide C (attacked)
                ('dgebf', other_epox_c, 'c1'),        # other opened epoxide C
            ],
            vect=(2.0, 0.0, 0.0),  # space the two reactants out a bit
            delete_byproduct=False,  # nothing leaves; H is conserved
        )
        return rxn, h_transfer, h_remaining, epox_o, other_epox_c

    def test_atom_count_unchanged(self):
        """Rxn 1 transfers an H — no atoms appear or disappear."""
        rxn, *_ = self._build()
        assert len(rxn.pre.mol.atoms) == len(rxn.post.mol.atoms)
        # And no byproduct cluster either.
        rxn.post.mol.clusters.compute_clusters()
        assert len(rxn.post.mol.clusters) == 1

    def test_h_transfer_topology(self):
        """The transferred H must move from N to O between pre and post."""
        rxn, h_transfer, h_remaining, epox_o, _ = self._build()
        n_new = rxn.id_map[('detda', DETDA_N_PRIMARY)]
        h_t_new = rxn.id_map[('detda', h_transfer)]
        h_r_new = rxn.id_map[('detda', h_remaining)]
        o_new = rxn.id_map[('dgebf', epox_o)]
        c_new = rxn.id_map[('dgebf', DGEBF_EPOXIDE_C)]

        # Pre: H_t bonded to N, NOT to O. Epoxide O bonded to attacked C.
        assert rxn.pre.mol.bonds.key_rule(n_new, h_t_new) in rxn.pre.mol.bonds
        assert rxn.pre.mol.bonds.key_rule(o_new, h_t_new) not in rxn.pre.mol.bonds
        assert rxn.pre.mol.bonds.key_rule(o_new, c_new) in rxn.pre.mol.bonds

        # Post: H_t bonded to O, NOT to N. Epoxide O-C bond gone.
        assert rxn.post.mol.bonds.key_rule(n_new, h_t_new) not in rxn.post.mol.bonds
        assert rxn.post.mol.bonds.key_rule(o_new, h_t_new) in rxn.post.mol.bonds
        assert rxn.post.mol.bonds.key_rule(o_new, c_new) not in rxn.post.mol.bonds

        # The OTHER amine H stays put on N in both states.
        assert rxn.pre.mol.bonds.key_rule(n_new, h_r_new) in rxn.pre.mol.bonds
        assert rxn.post.mol.bonds.key_rule(n_new, h_r_new) in rxn.post.mol.bonds

        # The new N-C covalent bond appears only in post.
        assert rxn.post.mol.bonds.key_rule(n_new, c_new) in rxn.post.mol.bonds
        assert rxn.pre.mol.bonds.key_rule(n_new, c_new) not in rxn.pre.mol.bonds

    def test_n_valence_preserved(self):
        """N had 3 bonds before (2H + 1 ring C); should still have 3 after
        (1H + 1 ring C + 1 new C)."""
        rxn, *_ = self._build()
        n_new = rxn.id_map[('detda', DETDA_N_PRIMARY)]
        pre_n_bonds = sum(1 for k in rxn.pre.mol.bonds if n_new in k)
        post_n_bonds = sum(1 for k in rxn.post.mol.bonds if n_new in k)
        assert pre_n_bonds == 3
        assert post_n_bonds == 3

    def test_o_valence_preserved(self):
        """Epoxide O had 2 bonds (to both epoxide C's); the alcohol O has 2
        bonds (to remaining epoxide C and to the transferred H)."""
        rxn, _, _, epox_o, _ = self._build()
        o_new = rxn.id_map[('dgebf', epox_o)]
        pre_o_bonds = sum(1 for k in rxn.pre.mol.bonds if o_new in k)
        post_o_bonds = sum(1 for k in rxn.post.mol.bonds if o_new in k)
        assert pre_o_bonds == 2
        assert post_o_bonds == 2

    def test_attacked_c_valence_preserved(self):
        """The attacked epoxide C had 4 bonds (O, other-C, methylene-C, H);
        post-reaction it should still have 4 (other-C, methylene-C, H, N)."""
        rxn, *_ = self._build()
        c_new = rxn.id_map[('dgebf', DGEBF_EPOXIDE_C)]
        pre_c_bonds = sum(1 for k in rxn.pre.mol.bonds if c_new in k)
        post_c_bonds = sum(1 for k in rxn.post.mol.bonds if c_new in k)
        assert pre_c_bonds == 4
        assert post_c_bonds == 4

    def test_nta_changes_applied(self):
        rxn, h_transfer, _, epox_o, other_epox_c = self._build()
        n_new = rxn.id_map[('detda', DETDA_N_PRIMARY)]
        h_new = rxn.id_map[('detda', h_transfer)]
        o_new = rxn.id_map[('dgebf', epox_o)]
        c_new = rxn.id_map[('dgebf', DGEBF_EPOXIDE_C)]
        c2_new = rxn.id_map[('dgebf', other_epox_c)]

        # N is intentionally NOT in changed_typ — keeps its 'nn' NTA.
        assert 'nn' in rxn.post.mol.atoms[n_new].comment
        assert rxn.post.mol.atoms[h_new].comment == 'ho'
        assert rxn.post.mol.atoms[o_new].comment == 'oh'
        assert rxn.post.mol.atoms[c_new].comment == 'c2'
        assert rxn.post.mol.atoms[c2_new].comment == 'c1'

        # Pre is untouched; original NTA strings still appear in the comment.
        assert 'nn' in rxn.pre.mol.atoms[n_new].comment
        assert 'hn2' in rxn.pre.mol.atoms[h_new].comment
        assert 'o3e' in rxn.pre.mol.atoms[o_new].comment
        assert 'c3h' in rxn.pre.mol.atoms[c_new].comment

    def test_remaining_h_stays_amine_h(self):
        """The amine H that does NOT transfer keeps its NTA in the post."""
        rxn, _, h_remaining, *_ = self._build()
        h_r_new = rxn.id_map[('detda', h_remaining)]
        assert 'hn2' in rxn.post.mol.atoms[h_r_new].comment

    def test_initiators_in_map(self):
        rxn, *_ = self._build()
        n_new = rxn.id_map[('detda', DETDA_N_PRIMARY)]
        c_new = rxn.id_map[('dgebf', DGEBF_EPOXIDE_C)]
        assert rxn.rxn_map['InitiatorIDs'] == [n_new, c_new]

    def test_initiator_separation_via_vect(self):
        rxn, *_ = self._build()
        n_new = rxn.id_map[('detda', DETDA_N_PRIMARY)]
        c_new = rxn.id_map[('dgebf', DGEBF_EPOXIDE_C)]
        n_atom = rxn.pre.mol.atoms[n_new]
        c_atom = rxn.pre.mol.atoms[c_new]
        # vect=(2.0, 0, 0): C placed at N + (2,0,0).
        assert c_atom.x == pytest.approx(n_atom.x + 2.0, abs=1e-9)
        assert c_atom.y == pytest.approx(n_atom.y, abs=1e-9)
        assert c_atom.z == pytest.approx(n_atom.z, abs=1e-9)


class TestReaction2SecondaryAmineEpoxide:
    """
    Rxn 2 — same mechanism on a secondary amine.

    We build the secondary-amine input by running Rxn 1 first and re-using its
    ``post`` Fragment. The N keeps its preserved id, so we can attack it again
    with another DGEBF.
    """

    def _build_post_rxn1(self):
        """Reuse Rxn 1 to produce a half-reacted post Fragment."""
        detda = _load_detda()
        dgebf = _load_dgebf()
        h_transfer = _amine_hs(detda, DETDA_N_PRIMARY)[0]
        epox_o, other_epox_c = _epoxide_partners(dgebf, DGEBF_EPOXIDE_C)

        A = Fragment.from_initiator_steps(detda, DETDA_N_PRIMARY, n_steps=DETDA_DEPTH)
        B = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)
        rxn1 = ReactionTemplate(
            A, B,
            A_name='detda', B_name='dgebf',
            breaks=[
                ('detda', DETDA_N_PRIMARY, h_transfer),
                ('dgebf', epox_o, DGEBF_EPOXIDE_C),
            ],
            makes=[
                ('detda', DETDA_N_PRIMARY, 'dgebf', DGEBF_EPOXIDE_C),
                ('detda', h_transfer, 'dgebf', epox_o),
            ],
            changed_typ=[
                ('detda', h_transfer, 'ho'),
                ('dgebf', epox_o, 'oh'),
                ('dgebf', DGEBF_EPOXIDE_C, 'c2'),
                ('dgebf', other_epox_c, 'c1'),
            ],
            vect=(2.0, 0.0, 0.0),
            delete_byproduct=False,
        )
        return rxn1.post, rxn1.id_map

    def test_secondary_amine_attacks_second_epoxide(self):
        # Rxn 1 product is the secondary-amine-bearing molecule.
        post1, id_map1 = self._build_post_rxn1()
        n_id = id_map1[('detda', DETDA_N_PRIMARY)]
        # The N is now secondary: 1 amine H, 1 ring C, 1 new C (from Rxn 1).
        graph_post1 = post1.mol.generate_graph()
        n_neighbors = graph_post1[n_id]
        amine_h_remaining = [nbr for nbr in n_neighbors
                             if 'hn2' in post1.mol.atoms[nbr].comment]
        assert len(amine_h_remaining) == 1, \
            f"expected exactly one amine H on secondary N, got {amine_h_remaining}"
        h_transfer_2 = amine_h_remaining[0]

        # Build a second DGEBF for the next epoxide attack.
        dgebf2 = _load_dgebf()
        epox_o2, other_epox_c2 = _epoxide_partners(dgebf2, DGEBF_EPOXIDE_C)

        A2 = Fragment.from_initiator_steps(post1, n_id, n_steps=DETDA_DEPTH)
        B2 = Fragment.from_initiator_steps(dgebf2, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)

        rxn2 = ReactionTemplate(
            A2, B2,
            A_name='detda_sec', B_name='dgebf2',
            breaks=[
                ('detda_sec', n_id, h_transfer_2),       # second N-H
                ('dgebf2', epox_o2, DGEBF_EPOXIDE_C),
            ],
            makes=[
                ('detda_sec', n_id, 'dgebf2', DGEBF_EPOXIDE_C),
                ('detda_sec', h_transfer_2, 'dgebf2', epox_o2),
            ],
            changed_typ=[
                # N stays 'nn' through both reactions in this FF practice.
                ('detda_sec', h_transfer_2, 'ho'),
                ('dgebf2', epox_o2, 'oh'),
                ('dgebf2', DGEBF_EPOXIDE_C, 'c2'),
                ('dgebf2', other_epox_c2, 'c1'),
            ],
            vect=(2.0, 0.0, 0.0),
            delete_byproduct=False,
        )

        # After Rxn 2, the N has: 1 ring C + 2 new C's (from Rxn 1 and Rxn 2) = 3 bonds, no H's.
        n_new = rxn2.id_map[('detda_sec', n_id)]
        post_n_neighbors = rxn2.post.mol.generate_graph()[n_new]
        h_neighbors = [nbr for nbr in post_n_neighbors
                       if 'hn2' in rxn2.post.mol.atoms[nbr].comment]
        assert len(h_neighbors) == 0, "tertiary amine should have no amine H"
        # N's NTA stays 'nn' (we don't tag tertiary amines separately here).
        assert 'nn' in rxn2.post.mol.atoms[n_new].comment


class TestValidation:
    def test_non_fragment_input_rejected(self):
        detda = _load_detda()
        with pytest.raises(TypeError):
            ReactionTemplate(detda, detda)

    def test_break_bond_not_found_raises(self):
        detda = _load_detda()
        dgebf = _load_dgebf()
        A = Fragment.from_initiator_steps(detda, DETDA_N_PRIMARY, n_steps=DETDA_DEPTH)
        B = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)
        with pytest.raises(ValueError):
            # Both 28 and 29 are amine Hs on N=12 but they are NOT bonded to
            # each other.
            ReactionTemplate(
                A, B,
                A_name='detda', B_name='dgebf',
                breaks=[('detda', 28, 29)],
            )


class TestFileIO:
    def test_round_trip_through_skeletal_data(self):
        detda = _load_detda()
        dgebf = _load_dgebf()
        amine_hs = _amine_hs(detda, DETDA_N_PRIMARY)
        epox_o, other_epox_c = _epoxide_partners(dgebf, DGEBF_EPOXIDE_C)
        A = Fragment.from_initiator_steps(detda, DETDA_N_PRIMARY, n_steps=DETDA_DEPTH)
        B = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)
        rxn = ReactionTemplate(
            A, B,
            A_name='detda', B_name='dgebf',
            breaks=[
                ('detda', DETDA_N_PRIMARY, amine_hs[0]),
                ('dgebf', epox_o, DGEBF_EPOXIDE_C),
            ],
            makes=[
                ('detda', DETDA_N_PRIMARY, 'dgebf', DGEBF_EPOXIDE_C),
                ('detda', amine_hs[0], 'dgebf', epox_o),
            ],
            delete_byproduct=False,
        )

        with tempfile.TemporaryDirectory() as td:
            pre = Path(td) / Path('pre.data')
            post = Path(td) / Path('post.data')
            map_path = Path(td) / Path('rxn.txt')
            rxn.write_data(pre, post)
            rxn.write_map(map_path)

            assert os.path.exists(pre) and os.path.exists(post) and os.path.exists(map_path)

            for path in (pre, post):
                roundtrip = Molspace(path)
                assert len(roundtrip.atoms) == len(rxn.pre.mol.atoms)
                assert len(roundtrip.angles) == 0
                assert len(roundtrip.dihedrals) == 0
                assert len(roundtrip.impropers) == 0
                for bond in roundtrip.bonds.values():
                    assert bond.type == 1

            assert 'mooonpy map file' in open(map_path).read()


class TestFragmentEdges:
    """Edge detection, single-interior-bond rule, severed-neighbor NTAs, and
    Molspace-or-Fragment input handling."""

    def test_boundary_equals_edges_and_single_interior_bond(self):
        dgebf = _load_dgebf()
        frag = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)
        graph = frag.mol.generate_graph()
        assert frag.edges
        for edge in frag.edges:
            n_interior = sum(1 for n in graph[edge] if n in frag.mol.atoms)
            assert n_interior == 1

    def test_edge_external_ntas_match_direct_graph_computation(self):
        detda = _load_detda()
        frag = Fragment.from_initiator_steps(detda, DETDA_N_PRIMARY, n_steps=DETDA_DEPTH)
        full_graph = detda.generate_graph()
        interior = set(frag.mol.atoms)
        for edge_id, ntas in frag.edge_external_ntas.items():
            externals = [detda.atoms[n].comment
                         for n in full_graph[edge_id] if n not in interior]
            assert sorted(ntas) == sorted(externals)
            assert ntas, f"edge {edge_id} should have severed neighbors"

    def test_from_interior_derives_boundary_edges(self):
        detda = _load_detda()
        interior = {DETDA_N_PRIMARY, 28, 29, 6}  # N, two amine H, ring C
        frag = Fragment.from_interior(detda, DETDA_N_PRIMARY, interior)
        graph = detda.generate_graph()
        expected = {i for i in interior
                    if any(n not in interior for n in graph[i])}
        assert frag.edges == expected
        # Amine H atoms have no external neighbors -> not edges.
        assert 28 not in frag.edges and 29 not in frag.edges

    def test_unreached_declared_edge_is_dropped(self):
        detda = _load_detda()
        frag = Fragment.from_initiator_and_edges(
            detda, initiator=DETDA_N_PRIMARY, edges=[6, 31], max_steps=2,
        )
        assert 6 in frag.edges
        assert 31 not in frag.edges
        assert 31 not in frag.edge_external_ntas

    def test_multi_interior_bond_edge_warns(self):
        """An epoxide-ring-only interior makes both ring C's edges with two
        interior bonds each; fix bond/react cannot anchor those."""
        dgebf = _load_dgebf()
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always')
            frag = Fragment.from_interior(dgebf, 2, {1, 2, 3})
            msgs = [str(w.message) for w in ws]
        assert any('bonds into the fragment interior' in m for m in msgs)
        assert frag.edges == {2, 3}

    def test_assembled_keeps_explicit_edge_ntas(self):
        dgebf = _load_dgebf()
        explicit = {6: ['cp', 'cp']}
        no_remap = {2, 3}
        frag = Fragment.assembled(dgebf, edges=[6],
                                  edge_external_ntas=explicit,
                                  no_remap_ids=no_remap)
        assert frag.initiator is None
        assert frag.edges == {6}
        assert frag.edge_external_ntas == explicit
        assert frag.no_remap_ids == no_remap
        assert len(frag.mol.atoms) == len(dgebf.atoms)

    def test_resolve_input_accepts_molspace_and_fragment(self):
        dgebf = _load_dgebf()
        from_mol = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=2)
        from_frag = Fragment.from_initiator_steps(from_mol, DGEBF_EPOXIDE_C, n_steps=2)
        assert set(from_frag.mol.atoms).issubset(set(from_mol.mol.atoms))

    def test_fragment_of_fragment_inherits_surviving_edge_ntas(self):
        dgebf = _load_dgebf()
        B = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)
        assert B.edges == {6}
        assert B.edge_external_ntas[6]  # non-empty

        # 6's external neighbors were pruned from B.mol, so the computed list
        # is empty and the constructor must fall back to B's stored info.
        B2 = Fragment.from_initiator_steps(B, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)
        assert 6 in B2.edges
        assert B2.edge_external_ntas[6] == B.edge_external_ntas[6]

    def test_fragment_of_fragment_new_edge_computed_fresh(self):
        dgebf = _load_dgebf()
        B = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)
        B2 = Fragment.from_initiator_steps(B, DGEBF_EPOXIDE_C, n_steps=3)
        new_edges = B2.edges - B.edges
        assert new_edges
        graph = B.mol.generate_graph()
        interior = set(B2.mol.atoms)
        for e in new_edges:
            expected = [B.mol.atoms[n].comment
                        for n in graph[e] if n not in interior]
            assert sorted(B2.edge_external_ntas[e]) == sorted(expected)


class TestReactionTemplateEdges:
    """Edge info survives the combine + renumber and lands in the .nta file."""

    def _build(self):
        detda = _load_detda()
        dgebf = _load_dgebf()
        h_transfer = _amine_hs(detda, DETDA_N_PRIMARY)[0]
        epox_o, other_epox_c = _epoxide_partners(dgebf, DGEBF_EPOXIDE_C)
        A = Fragment.from_initiator_steps(detda, DETDA_N_PRIMARY, n_steps=DETDA_DEPTH)
        B = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)
        rxn = ReactionTemplate(
            A, B,
            A_name='detda', B_name='dgebf',
            breaks=[
                ('detda', DETDA_N_PRIMARY, h_transfer),
                ('dgebf', epox_o, DGEBF_EPOXIDE_C),
            ],
            makes=[
                ('detda', DETDA_N_PRIMARY, 'dgebf', DGEBF_EPOXIDE_C),
                ('detda', h_transfer, 'dgebf', epox_o),
            ],
            changed_typ=[
                ('detda', h_transfer, 'ho'),
                ('dgebf', epox_o, 'oh'),
                ('dgebf', DGEBF_EPOXIDE_C, 'c2'),
                ('dgebf', other_epox_c, 'c1'),
            ],
            delete_byproduct=False,
        )
        return rxn, A, B

    def test_pre_and_post_are_fragments(self):
        rxn, _, _ = self._build()
        assert isinstance(rxn.pre, Fragment)
        assert isinstance(rxn.post, Fragment)
        assert rxn.pre.initiator is None  # assembled whole-template form

    def test_edge_ids_translated_to_renumbered_space(self):
        rxn, A, B = self._build()
        expected = sorted(
            [rxn.id_map[('detda', e)] for e in A.edges] +
            [rxn.id_map[('dgebf', e)] for e in B.edges]
        )
        assert rxn.rxn_map['EdgeIDs'] == expected
        for eid in rxn.rxn_map['EdgeIDs']:
            assert eid in rxn.pre.mol.atoms

    def test_pre_post_edge_external_ntas_identical(self):
        rxn, _, _ = self._build()
        assert rxn.pre.edge_external_ntas == rxn.post.edge_external_ntas

    def test_edge_external_ntas_keyed_in_renumbered_space(self):
        rxn, A, B = self._build()
        combined = {}
        for e, ntas in A.edge_external_ntas.items():
            combined[rxn.id_map[('detda', e)]] = sorted(ntas)
        for e, ntas in B.edge_external_ntas.items():
            combined[rxn.id_map[('dgebf', e)]] = sorted(ntas)
        got = {k: sorted(v) for k, v in rxn.pre.edge_external_ntas.items()}
        assert got == combined

    def test_nta_file_has_edge_section(self):
        rxn, _, _ = self._build()
        with tempfile.TemporaryDirectory() as td:
            pre = Path(td) / Path('pre.data')
            post = Path(td) / Path('post.data')
            _, _, pre_nta, _ = rxn.write_data(pre, post)
            txt = open(pre_nta).read()
            assert 'edge id' in txt
            section = txt.split('edge id', 1)[1].strip().splitlines()
            assert section
            for line in section:
                parts = line.split()
                eid = int(parts[0])
                assert eid in rxn.rxn_map['EdgeIDs']
                assert len(parts) > 1

    def test_post_fragment_feeds_next_reaction(self):
        rxn, _, _ = self._build()
        n_new = rxn.id_map[('detda', DETDA_N_PRIMARY)]
        nxt = Fragment.from_initiator_steps(rxn.post, n_new, n_steps=2)
        assert isinstance(nxt, Fragment)
        assert set(nxt.mol.atoms).issubset(set(rxn.post.mol.atoms))

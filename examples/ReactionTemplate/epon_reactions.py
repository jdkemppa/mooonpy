# -*- coding: utf-8 -*-
"""
Build the two EPON 862 (DETDA + DGEBF) crosslink reaction templates and write
them to the same folder as this script.

Reaction 1 — primary amine + epoxide -> secondary amine + alcohol
    R-N(H)(H) + C(H)(R')-O-C(H)(R'')   ->   R-N(H)-C(H)(R')-C(OH)(H)(R'')

Reaction 2 — secondary amine + epoxide -> tertiary amine + alcohol
    Same mechanism, run on the post-Rxn1 molecule (the N now has only one H).

For both reactions the amine N attacks one epoxide C, the C-O bond on the
attacked side breaks, and one amine H transfers from N to the now-alcohol O.
Atom ids are preserved 1:1 between pre and post.

Outputs (eight files, plus this script's own logs to stdout):
    pre_rxn1.data, post_rxn1.data, pre_rxn1.nta, post_rxn1.nta, rxn1_map.txt
    pre_rxn2.data, post_rxn2.data, pre_rxn2.nta, post_rxn2.nta, rxn2_map.txt

Run with ``python -i epon_reactions.py`` to keep all module-level variables
(``detda``, ``dgebf``, ``A``, ``B``, ``rxn1``, ``rxn2``, ...) live for
inspection afterward.
"""
import os
import sys

# Allow running this script in-place against the dev worktree without an
# editable install. Remove after merge if mooonpy is on PYTHONPATH globally.
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', 'src')
))

import mooonpy
from mooonpy import Molspace, Path, ReactionTemplate
from mooonpy.template import Fragment, canonicalize_ntas


HERE = Path(os.path.dirname(os.path.abspath(__file__)))
EPON_DIR = HERE / Path('../EPON_862/atom_typing_Outputs')
OUT_DIR = HERE  # write outputs alongside this script


# ----------------------------------------------------------------------------
# Reactive-site IDs (by inspection of the modern atom_typing-output monomers).
# ----------------------------------------------------------------------------
DETDA_N_PRIMARY = 12       # one of two primary amines (the other is N=13)
DGEBF_EPOXIDE_C = 2        # epoxide C the amine attacks

# Fragment depth: edges must sit ≥3 bonds out from the reactive atom so every
# dihedral that crosses the new/broken bond has all four atoms inside the
# template (4 atoms = 3 bonds), letting fix bond/react match parameters.
# DGEBF gets one extra step to comfortably reach the phenyl ring on the other
# side of the methylene-O linker.
DETDA_DEPTH = 3
DGEBF_DEPTH = 4


def amine_hs(mol, n_id):
    """Walk one bond from N to find the amine H atoms (NTA = 'hn2')."""
    return [nbr for nbr in mol.generate_graph()[n_id]
            if mol.atoms[nbr].comment == 'hn2']


def epoxide_partners(mol, attacked_c_id):
    """From the attacked epoxide C, find the O and the other epoxide C."""
    nbrs = mol.generate_graph()[attacked_c_id]
    o_id = next(nbr for nbr in nbrs if mol.atoms[nbr].comment == 'o3e')
    other_c = next(nbr for nbr in nbrs
                   if mol.atoms[nbr].comment == 'c3h' and nbr != attacked_c_id)
    return o_id, other_c


def amine_h_on(mol, n_id):
    """Secondary amine N has exactly one remaining 'hn2' neighbor."""
    hs = amine_hs(mol, n_id)
    if len(hs) != 1:
        raise RuntimeError(
            f"expected exactly one amine H on N={n_id}, found {hs}"
        )
    return hs[0]


if __name__ == '__main__':
    print(f'mooonpy from: {mooonpy.__file__}')
    print(f'writing outputs to: {OUT_DIR}')
    print()

    # ------------------------------------------------------------------
    # Load monomers and canonicalize NTAs first.
    # atom_typing writes verbose comments like
    #   "C        Sp2     avg-angle: 120.000000      type: cp"
    # canonicalize_ntas() rewrites every atom.comment to just the NTA token
    # ('cp'), so downstream lookups, .data writes, and .nta writes all see
    # the same clean strings.
    # ------------------------------------------------------------------
    detda = Molspace(EPON_DIR / Path('detda_typed.data'))
    dgebf = Molspace(EPON_DIR / Path('dgebf_typed.data'))
    canonicalize_ntas(detda)
    canonicalize_ntas(dgebf)
    print(f'loaded DETDA: {len(detda.atoms)} atoms, masses={list(detda.ff.masses.keys())}')
    print(f'loaded DGEBF: {len(dgebf.atoms)} atoms, masses={list(dgebf.ff.masses.keys())}')
    print()

    # ------------------------------------------------------------------
    # Reaction 1 — primary amine + epoxide
    # ------------------------------------------------------------------
    print('=== Reaction 1: primary amine + epoxide ===')
    h_transfer = amine_hs(detda, DETDA_N_PRIMARY)[0]
    epox_o, other_epox_c = epoxide_partners(dgebf, DGEBF_EPOXIDE_C)
    print(f'  N initiator        : {DETDA_N_PRIMARY}')
    print(f'  amine H to transfer: {h_transfer}')
    print(f'  epoxide O / other C: {epox_o} / {other_epox_c}')

    A = Fragment.from_initiator_steps(detda, DETDA_N_PRIMARY, n_steps=DETDA_DEPTH)
    B = Fragment.from_initiator_steps(dgebf, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)

    rxn1 = ReactionTemplate(
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
        # N keeps NTA 'nn' through both reactions (we don't tag primary,
        # secondary, tertiary amines separately in this FF practice).
        changed_typ=[
            ('detda', h_transfer, 'ho'),          # amine H -> hydroxyl H
            ('dgebf', epox_o, 'oh'),              # epoxide O -> hydroxyl O
            ('dgebf', DGEBF_EPOXIDE_C, 'c2'),     # opened epoxide C (the attacked one)
            ('dgebf', other_epox_c, 'c1'),        # other opened epoxide C
        ],
        vect=(2.0, 3.0, 1.0),
        delete_byproduct=False,                   # H is conserved; no byproduct
    )

    pre1 = OUT_DIR / Path('pre_rxn1.data')
    post1 = OUT_DIR / Path('post_rxn1.data')
    map1 = OUT_DIR / Path('rxn1_map.txt')
    pre1, post1, pre1_nta, post1_nta = rxn1.write_data(pre1, post1)
    rxn1.write_map(map1)

    print(f'  pre/post atom count: {len(rxn1.pre.mol.atoms)} / {len(rxn1.post.mol.atoms)}')
    print(f'  combined masses    : {list(rxn1.pre.mol.ff.masses.keys())}')
    print(f'  edge external NTAs : {rxn1.pre.edge_external_ntas}')
    print(f'  initiator ids      : {rxn1.rxn_map["InitiatorIDs"]}')
    print(f'  edge ids           : {rxn1.rxn_map["EdgeIDs"]}')
    print(f'  wrote              : {pre1.basename()}, {post1.basename()},')
    print(f'                       {pre1_nta.basename()}, {post1_nta.basename()},')
    print(f'                       {map1.basename()}')
    print()

    # ------------------------------------------------------------------
    # Reaction 2 — secondary amine + epoxide (uses Rxn 1's product as input)
    # ------------------------------------------------------------------
    print('=== Reaction 2: secondary amine + epoxide ===')
    n_in_post1 = rxn1.id_map[('detda', DETDA_N_PRIMARY)]
    h_transfer_2 = amine_h_on(rxn1.post.mol, n_in_post1)
    print(f'  N initiator (preserved id): {n_in_post1}')
    print(f'  remaining amine H         : {h_transfer_2}')

    dgebf2 = Molspace(EPON_DIR / Path('dgebf_typed.data'))
    canonicalize_ntas(dgebf2)
    epox_o2, other_epox_c2 = epoxide_partners(dgebf2, DGEBF_EPOXIDE_C)

    # rxn1.post is a Fragment; from_initiator_steps inherits its
    # edge_external_ntas for any edge that survives into this new slice.
    A2 = Fragment.from_initiator_steps(rxn1.post, n_in_post1, n_steps=5)
    B2 = Fragment.from_initiator_steps(dgebf2, DGEBF_EPOXIDE_C, n_steps=DGEBF_DEPTH)
    B2.mol.atoms.move((-1, 1, 1), mode='scale')
    rxn2 = ReactionTemplate(
        A2, B2,
        A_name='detda_sec', B_name='dgebf2',
        breaks=[
            ('detda_sec', n_in_post1, h_transfer_2),
            ('dgebf2', epox_o2, DGEBF_EPOXIDE_C),
        ],
        makes=[
            ('detda_sec', n_in_post1, 'dgebf2', DGEBF_EPOXIDE_C),
            ('detda_sec', h_transfer_2, 'dgebf2', epox_o2),
        ],
        changed_typ=[
            ('detda_sec', h_transfer_2, 'ho'),
            ('dgebf2', epox_o2, 'oh'),
            ('dgebf2', DGEBF_EPOXIDE_C, 'c2'),
            ('dgebf2', other_epox_c2, 'c1'),
        ],
        vect=(2.0, 0.0, -3.0),
        delete_byproduct=False,
    )

    pre2 = OUT_DIR / Path('pre_rxn2.data')
    post2 = OUT_DIR / Path('post_rxn2.data')
    map2 = OUT_DIR / Path('rxn2_map.txt')
    pre2, post2, pre2_nta, post2_nta = rxn2.write_data(pre2, post2)
    rxn2.write_map(map2)

    print(f'  pre/post atom count: {len(rxn2.pre.mol.atoms)} / {len(rxn2.post.mol.atoms)}')
    print(f'  combined masses    : {list(rxn2.pre.mol.ff.masses.keys())}')
    print(f'  edge external NTAs : {rxn2.pre.edge_external_ntas}')
    print(f'  initiator ids      : {rxn2.rxn_map["InitiatorIDs"]}')
    print(f'  edge ids           : {rxn2.rxn_map["EdgeIDs"]}')
    print(f'  wrote              : {pre2.basename()}, {post2.basename()},')
    print(f'                       {pre2_nta.basename()}, {post2_nta.basename()},')
    print(f'                       {map2.basename()}')
    print()

    # ------------------------------------------------------------------
    # Sanity: confirm the H really moved from N to O in Rxn 1's post.
    # ------------------------------------------------------------------
    print('=== Sanity: H-transfer topology (Rxn 1) ===')
    h_new = rxn1.id_map[('detda', h_transfer)]
    pre_h_bonds = [k for k in rxn1.pre.mol.bonds if h_new in k]
    post_h_bonds = [k for k in rxn1.post.mol.bonds if h_new in k]
    print(f'  H atom id in template       : {h_new}')
    print(f'  bonds involving H in pre    : {pre_h_bonds}')
    print(f'  bonds involving H in post   : {post_h_bonds}')
    print(f'  H NTA in post               : {rxn1.post.mol.atoms[h_new].comment!r}')

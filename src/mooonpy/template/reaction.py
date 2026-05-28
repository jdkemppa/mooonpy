# -*- coding: utf-8 -*-
"""
ReactionTemplate — build a LAMMPS ``fix bond/react`` template pair from two
:class:`Fragment` objects.

Atom ids are preserved 1:1 between the pre and post Molspaces so the LAMMPS
internal topology-map optimizer is sidestepped.

Force-field handling stays out of this module: the output ``.data`` files are
emitted with no angles / dihedrals / impropers, all bonds collapsed to a single
dummy bond type, and FF coefficient dicts cleared. NTA strings live in
``atom.comment`` and are the truth source for downstream typing tools, which
overwrite the integer ``atom.type`` regardless.
"""
import warnings

from mooonpy.tools.file_utils import Path
from mooonpy.template.fragment import Fragment
from mooonpy.template.mapfile import write_mapfile as _write_mapfile_fn
from mooonpy.template.ntafile import (
    write_ntafile as _write_ntafile_fn,
    canonicalize_ntas as _canonicalize_ntas,
)


class ReactionTemplate:
    """
    Build a ``fix bond/react`` template pair from two :class:`Fragment` reactants.

    :param A: First reactant fragment.
    :type A: Fragment
    :param B: Second reactant fragment.
    :type B: Fragment
    :param A_name: Short string used as the dictionary key for ``A`` atoms in
                   ``breaks``, ``makes``, ``changed_typ``, and ``id_map``.
    :type A_name: str
    :param B_name: Same as ``A_name`` for ``B``.
    :type B_name: str
    :param breaks: List of ``(name, orig_id1, orig_id2)`` triples — bonds to
                   remove from the post state. Both ids must reference the
                   same reactant.
    :type breaks: list[tuple[str, int, int]]
    :param makes: List of ``(name1, orig_id1, name2, orig_id2)`` quads — bonds
                  to create in the post state.
    :type makes: list[tuple[str, int, str, int]]
    :param changed_typ: List of ``(name, orig_id, new_nta)`` triples — write
                        ``new_nta`` into ``post.atoms[new_id].comment``. The
                        integer ``atom.type`` is left untouched (downstream
                        atom_typing re-assigns it from the ``.nta`` file).
    :type changed_typ: list[tuple[str, int, str]]
    :param vect: Relative displacement between the two fragment initiators in
                 the combined pre state. ``B``'s initiator is placed at
                 ``A_initiator + vect``. ``(0, 0, 0)`` makes them coincident.
    :type vect: tuple[float, float, float]
    :param delete_byproduct: After applying bond changes, run cluster ID on
                             the post state; everything not in the largest
                             cluster is added to ``DeleteIDs``.
    :type delete_byproduct: bool

    Computed attributes:

    * ``pre`` (Molspace): merged pre-reaction state with 1..N ids.
    * ``post`` (Molspace): same ids; bonds and NTAs updated.
    * ``rxn_map`` (dict): editable intermediate; written by :meth:`write_map`.
    * ``id_map`` (dict): ``{(name, original_id): preserved_id}``.
    """

    def __init__(self, A, B, *,
                 A_name='A', B_name='B',
                 breaks=(), makes=(), changed_typ=(),
                 vect=(0.0, 0.0, 0.0),
                 delete_byproduct=True):
        self._validate_inputs(A, B, A_name, B_name, breaks, makes, changed_typ)

        self.A_name = A_name
        self.B_name = B_name
        self.A = A
        self.B = B

        # Work on copies of the fragment Molspaces so we don't mutate inputs.
        A_mol = A.mol.copy()
        B_mol = B.mol.copy()

        # Place B's initiator at A's initiator + vect (relative-offset positioning).
        a_init = A.initiator
        b_init = B.initiator
        a_atom = A_mol.atoms[a_init]
        b_atom = B_mol.atoms[b_init]
        dx = a_atom.x + float(vect[0]) - b_atom.x
        dy = a_atom.y + float(vect[1]) - b_atom.y
        dz = a_atom.z + float(vect[2]) - b_atom.z
        B_mol.atoms.move((dx, dy, dz), mode='offset')

        # Combine into a single pre state. shift_atomid is the offset applied
        # to B's atom ids during the combine. offset_types=True forces the
        # type tables to be merged side-by-side; condense_element below then
        # collapses duplicates per element. shrinkwrap fits the box around
        # the merged geometry with a small pad so neither all2lmp nor LAMMPS
        # complains about atoms outside the box.
        shift_atomid = A_mol.combine(B_mol, vect=None, offset_types=True,
                                     assign_clusters=False)
        if A_mol.ff.masses:
            A_mol.ff.condense_element(A_mol.atoms)
        A_mol.atoms.shrinkwrap(pad=5.0)

        # Renumber the combined system 1..N and build the id_map.
        old2new = A_mol.contiguous()
        id_map = {}
        for orig_id in A.mol.atoms:
            id_map[(A_name, orig_id)] = old2new[orig_id]
        for orig_id in B.mol.atoms:
            id_map[(B_name, orig_id)] = old2new[orig_id + shift_atomid]
        self.id_map = id_map

        pre_mol = A_mol
        post_mol = pre_mol.copy()

        # Apply breaks.
        for entry in breaks:
            name, oid1, oid2 = entry
            nid1 = id_map[(name, oid1)]
            nid2 = id_map[(name, oid2)]
            key = post_mol.bonds.key_rule(nid1, nid2)
            if key not in post_mol.bonds:
                raise ValueError(
                    f"break bond ({name} {oid1}, {oid2}) -> ({nid1}, {nid2}) "
                    f"not found in post.bonds"
                )
            del post_mol.bonds[key]

        # Apply makes.
        for entry in makes:
            name1, oid1, name2, oid2 = entry
            nid1 = id_map[(name1, oid1)]
            nid2 = id_map[(name2, oid2)]
            key = post_mol.bonds.key_rule(nid1, nid2)
            if key in post_mol.bonds:
                raise ValueError(
                    f"make bond ({name1} {oid1}, {name2} {oid2}) -> "
                    f"({nid1}, {nid2}) already exists in post.bonds"
                )
            bond = post_mol.bonds.bond_factory()
            bond.ordered = [nid1, nid2]
            bond.type = 1  # FF-stripped output uses a single dummy bond type
            post_mol.bonds.id_set(bond, nid1, nid2)

        # Apply changed_typ: overwrite atom.comment with the new NTA string.
        for entry in changed_typ:
            name, oid, new_nta = entry
            nid = id_map[(name, oid)]
            post_mol.atoms[nid].comment = new_nta

        pre_ntas = {a.comment for a in pre_mol.atoms.values()}
        post_ntas = {a.comment for a in post_mol.atoms.values()}
        # if pre_ntas != post_ntas:
        #     warnings.warn(
        #         "ReactionTemplate: pre and post have different NTA sets; "
        #         "downstream type-merging tools will need to reconcile types "
        #         "across all templates and monomers.",
        #         stacklevel=2,
        #     )

        # Translate per-fragment edge sets, severed-neighbor NTAs, and
        # no-remap sets into the renumbered space. The external-neighbor NTAs
        # describe atoms *outside* the template, so they are identical for
        # pre and post (changed_typ only touches in-template atoms).
        edge_ids = sorted(
            [id_map[(A_name, oid)] for oid in A.edges] +
            [id_map[(B_name, oid)] for oid in B.edges]
        )
        edge_external_ntas = {}
        for oid, ntas in A.edge_external_ntas.items():
            edge_external_ntas[id_map[(A_name, oid)]] = list(ntas)
        for oid, ntas in B.edge_external_ntas.items():
            edge_external_ntas[id_map[(B_name, oid)]] = list(ntas)

        no_remap = (
            {id_map[(A_name, oid)] for oid in A.no_remap_ids} |
            {id_map[(B_name, oid)] for oid in B.no_remap_ids}
        )
        # qedge: atoms whose charges *do* get remapped — everything not marked
        # to be preserved by either fragment.
        if no_remap:
            qedge = sorted(aid for aid in post_mol.atoms if aid not in no_remap)
        else:
            qedge = []

        # Byproduct detection.
        delete_ids = []
        if delete_byproduct:
            post_mol.clusters.compute_clusters()
            if len(post_mol.clusters) > 1:
                for molid, cluster in post_mol.clusters.items():
                    if molid == 1:
                        continue  # largest cluster (per Clusters.compute_clusters sort)
                    delete_ids.extend(sorted(cluster.members))

        initiator_new = [id_map[(A_name, A.initiator)],
                         id_map[(B_name, B.initiator)]]

        # pre and post are themselves Fragments so the edge / severed-neighbor
        # info travels with them — e.g. for write_data's .nta edge section, or
        # when a follow-up reaction slices the post for the next crosslink.
        self.pre = Fragment.assembled(
            pre_mol, edges=edge_ids,
            edge_external_ntas=edge_external_ntas, no_remap_ids=no_remap,
        )
        self.post = Fragment.assembled(
            post_mol, edges=edge_ids,
            edge_external_ntas=edge_external_ntas, no_remap_ids=no_remap,
        )

        self.rxn_map = {
            'Comment': f'nicknames {A_name} {B_name}',
            'InitiatorIDs': initiator_new,
            'EdgeIDs': edge_ids,
            'DeleteIDs': sorted(delete_ids),
            'Equivalences': [(aid, aid) for aid in sorted(post_mol.atoms)],
            'CustomChargesQedge': qedge,
        }

    # ---- IO ----

    def write_map(self, path):
        """
        Write ``self.rxn_map`` to ``path`` in mooonpy map-file format.

        :param path: Output path. Wrapped in :class:`Path` if it isn't one.
        :type path: str or :class:`mooonpy.tools.file_utils.Path`
        """
        path = Path(path)
        _write_mapfile_fn(path, self.rxn_map)

    def write_data(self, pre_path, post_path):
        """
        Write pre and post ``.data`` files plus their companion ``.nta`` files.

        Data files are FF-stripped: no angles / dihedrals / impropers, no FF
        coefficient sections, no Velocities section, and a single dummy bond
        type. The Masses section is kept so downstream tools can sanity-check
        atom types. Atom comments are reduced to just the NTA token.

        The ``.nta`` files are derived from the data file paths by replacing
        the extension. They feed ``atom_typing``/``all2lmp`` downstream.

        :param pre_path: Output path for the pre-reaction ``.data`` file.
        :param post_path: Output path for the post-reaction ``.data`` file.
        :type pre_path: str or :class:`Path`
        :type post_path: str or :class:`Path`
        :return: ``(pre_data, post_data, pre_nta, post_nta)`` as :class:`Path` objects.
        """
        pre_path = Path(pre_path)
        post_path = Path(post_path)
        pre_nta = pre_path.new_ext('.nta')
        post_nta = post_path.new_ext('.nta')
        _write_skeletal_data(self.pre.mol, pre_path)
        _write_skeletal_data(self.post.mol, post_path)
        _write_ntafile_fn(pre_nta, self.pre.mol, self.pre.edge_external_ntas)
        _write_ntafile_fn(post_nta, self.post.mol, self.post.edge_external_ntas)
        return pre_path, post_path, pre_nta, post_nta

    # ---- validation ----

    @staticmethod
    def _validate_inputs(A, B, A_name, B_name, breaks, makes, changed_typ):
        if not isinstance(A, Fragment):
            raise TypeError("A must be a Fragment")
        if not isinstance(B, Fragment):
            raise TypeError("B must be a Fragment")
        if not isinstance(A_name, str):
            raise TypeError("A_name must be str")
        if not isinstance(B_name, str):
            raise TypeError("B_name must be str")
        if A_name == B_name:
            raise ValueError("A_name and B_name must differ")

        for entry in breaks:
            if len(entry) != 3:
                raise ValueError(f"break entry must be (name, id, id), got {entry!r}")
            name, oid1, oid2 = entry
            if name not in (A_name, B_name):
                raise ValueError(f"break name {name!r} is not A_name or B_name")
            if not isinstance(oid1, int) or not isinstance(oid2, int):
                raise TypeError(f"break ids must be int, got {entry!r}")
        for entry in makes:
            if len(entry) != 4:
                raise ValueError(f"make entry must be (name, id, name, id), got {entry!r}")
            name1, oid1, name2, oid2 = entry
            for nm in (name1, name2):
                if nm not in (A_name, B_name):
                    raise ValueError(f"make name {nm!r} is not A_name or B_name")
            if not isinstance(oid1, int) or not isinstance(oid2, int):
                raise TypeError(f"make ids must be int, got {entry!r}")
        for entry in changed_typ:
            if len(entry) != 3:
                raise ValueError(f"changed_typ entry must be (name, id, nta), got {entry!r}")
            name, oid, nta = entry
            if name not in (A_name, B_name):
                raise ValueError(f"changed_typ name {name!r} is not A_name or B_name")
            if not isinstance(oid, int):
                raise TypeError(f"changed_typ id must be int, got {entry!r}")
            if not isinstance(nta, str):
                raise TypeError(f"changed_typ nta must be str, got {entry!r}")


def _write_skeletal_data(mol, path):
    """
    Write ``mol`` as a LAMMPS .data file with all force-field info stripped.

    On a temporary copy:
      * Angles/dihedrals/impropers are dropped.
      * Every bond is collapsed to type 1.
      * All FF coefficient dicts are cleared (Masses survives).
      * Velocities are zeroed (the writer suppresses an all-zero Velocities
        section automatically).
      * Atom comments are canonicalized to just the NTA token, so downstream
        diffs and the companion .nta file stay clean.

    The original Molspace is not mutated.
    """
    skel = mol.copy()
    skel.angles.clear()
    skel.dihedrals.clear()
    skel.impropers.clear()
    for bond in skel.bonds.values():
        bond.type = 1
    skel.ff.pair_coeffs.clear()
    skel.ff.bond_coeffs.clear()
    skel.ff.angle_coeffs.clear()
    skel.ff.dihedral_coeffs.clear()
    skel.ff.improper_coeffs.clear()
    skel.ff.bondbond_coeffs.clear()
    skel.ff.bondangle_coeffs.clear()
    skel.ff.angleangletorsion_coeffs.clear()
    skel.ff.endbondtorsion_coeffs.clear()
    skel.ff.middlebondtorsion_coeffs.clear()
    skel.ff.bondbond13_coeffs.clear()
    skel.ff.angletorsion_coeffs.clear()
    skel.ff.angleangle_coeffs.clear()
    skel.ff.has_type_labels = False
    for atom in skel.atoms.values():
        atom.vx = 0.0
        atom.vy = 0.0
        atom.vz = 0.0
    _canonicalize_ntas(skel)
    skel.write_files(path, atom_style='full')

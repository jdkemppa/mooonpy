# -*- coding: utf-8 -*-
"""
Atom-ID remapping primitives for Molspace.

A single primitive ``remap_ids(mol, old2new)`` handles three use cases:

1. Pure renumbering: when ``old2new`` covers every atom id in ``mol``.
2. Slicing: when ``old2new`` covers only a subset; the missing ids and any
   bond/angle/dihedral/improper entry referencing them are dropped.
3. Combined slice + renumber: any subset, mapped to any new ids.

The Python identities of ``mol.atoms``, ``mol.bonds``, ``mol.angles``,
``mol.dihedrals``, ``mol.impropers``, and each surviving ``Atom``/``Bond``/...
object are preserved. Only the dict keys (and ``obj.id`` / ``obj.ordered``)
change.

The helper ``build_old2new`` constructs the mapping for the common cases
(contiguous renumber, fixed offset, identity-with-keep-set).

Bond-graph reachability lives in :mod:`mooonpy.molspace.graph_theory.interface`
(``find_reachable_neighbors``, ``find_reachable_neighbors_all``,
``steps_from_initiator``).
"""


def _rebuild_atoms(atoms_dict, old2new):
    """
    Rebuild ``atoms_dict`` in place: remap keys, update ``atom.id``, drop
    atoms missing from ``old2new``. Atom object identities are preserved.
    """
    snapshot = list(atoms_dict.items())
    atoms_dict.clear()
    for old_id, atom in snapshot:
        new_id = old2new.get(old_id)
        if new_id is None:
            continue
        atom.id = new_id
        atoms_dict[new_id] = atom
    # Re-insert sorted by new id for deterministic iteration order
    items = sorted(atoms_dict.items())
    atoms_dict.clear()
    atoms_dict.update(items)


def _rebuild_topology(container, old2new):
    """
    Rebuild a Bonds/Angles/Dihedrals/Impropers container in place: remap each
    entry's ``ordered`` ids via ``old2new``, drop entries that reference any
    id missing from the map, re-insert using the container's own
    ``key_rule`` via ``id_set``. Topology-object identities are preserved.
    """
    snapshot = list(container.items())
    container.clear()
    for _old_key, obj in snapshot:
        if any(i_ not in old2new for i_ in obj.ordered):
            continue  # slice: any missing atom drops this entry
        new_ordered = [old2new[i_] for i_ in obj.ordered]
        obj.ordered = new_ordered
        container.id_set(obj, *new_ordered)
    items = sorted(container.items())
    container.clear()
    container.update(items)


def remap_ids(mol, old2new):
    """
    Renumber and/or slice ``mol`` in place.

    Parameters
    ----------
    mol : Molspace
    old2new : dict[int, int]
        Map from current atom id to new atom id. Any current atom whose id is
        not in this dict is dropped, and any topology entry referencing a
        dropped id is also dropped.

    Notes
    -----
    Container Python identities (``mol.atoms``, ``mol.bonds``, etc.) and
    surviving per-atom/per-bond object identities are preserved. Clusters are
    recomputed from scratch.
    """
    _rebuild_atoms(mol.atoms, old2new)
    _rebuild_topology(mol.bonds, old2new)
    _rebuild_topology(mol.angles, old2new)
    _rebuild_topology(mol.dihedrals, old2new)
    _rebuild_topology(mol.impropers, old2new)
    mol.clusters.compute_clusters()


def build_old2new(mol, mode='contiguous', offset=0, start=1, keep=None):
    """
    Build an ``old2new`` mapping (and its inverse) for use with ``remap_ids``.

    Parameters
    ----------
    mol : Molspace
    mode : {'contiguous', 'offset', 'identity'}
        - ``'contiguous'``: assign ``start, start+1, ...`` in current atom
          insertion order. With ``keep`` set, only the kept ids are renumbered
          (so calling ``remap_ids`` with this map slices to those ids).
        - ``'offset'``: ``new = old + offset`` for every atom (or each id in
          ``keep`` if given).
        - ``'identity'``: each id maps to itself. Most useful with ``keep`` to
          slice without renumbering.
    offset : int
        Used by ``'offset'`` mode.
    start : int
        Used by ``'contiguous'`` mode.
    keep : iterable[int] or None
        If given, restricts the mapping to these atom ids. Missing ids will be
        dropped if this map is passed to ``remap_ids``.

    Returns
    -------
    (old2new, new2old) : tuple[dict[int, int], dict[int, int]]
    """
    keep_set = set(keep) if keep is not None else None

    if mode == 'contiguous':
        old2new = {}
        new = start
        for old in mol.atoms:  # preserves insertion order
            if keep_set is not None and old not in keep_set:
                continue
            old2new[old] = new
            new += 1
    elif mode == 'offset':
        if keep_set is not None:
            old2new = {old: old + offset for old in mol.atoms if old in keep_set}
        else:
            old2new = {old: old + offset for old in mol.atoms}
    elif mode == 'identity':
        if keep_set is not None:
            old2new = {old: old for old in mol.atoms if old in keep_set}
        else:
            old2new = {old: old for old in mol.atoms}
    else:
        raise ValueError(
            f"Unknown mode {mode!r}; expected 'contiguous', 'offset', or 'identity'"
        )

    new2old = {v: k for k, v in old2new.items()}
    return old2new, new2old


def contiguous(mol, start=1):
    """
    Renumber atoms 1..N (or start..start+N-1) in current insertion order.

    Returns
    -------
    old2new : dict[int, int]
    """
    old2new, _ = build_old2new(mol, mode='contiguous', start=start)
    remap_ids(mol, old2new)
    return old2new


def slice_mol(mol, keep_ids, renumber=True):
    """
    In-place slice ``mol`` to retain only ``keep_ids``.

    Parameters
    ----------
    keep_ids : iterable[int]
    renumber : bool
        If True, surviving atoms are renumbered contiguously from 1. If False,
        original ids are kept (`'identity'` mode).

    Returns
    -------
    old2new : dict[int, int]
    """
    mode = 'contiguous' if renumber else 'identity'
    old2new, _ = build_old2new(mol, mode=mode, keep=keep_ids)
    remap_ids(mol, old2new)
    return old2new



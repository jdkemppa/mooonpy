# -*- coding: utf-8 -*-
import mooonpy.molspace.rings as rings
import mooonpy.molspace.force_field as FF
from mooonpy.molspace.graph_theory.ring_analysis import _reduced_graph

"""
Very rough first iteration of this code, I'll first assume that it coarse grains
by ring centers, but the eventual goal will be to apply some sort of reacter-like
template that is then converted into a particle at that template's center of mass.

...that requires built-in atom_typing-like functionality and an expanded graph theory library
"""


# depth first search algorighm to find connected rings within the domain of rgraph
def _dfs_chain(rgraph, graph, centers, start):
    """
    rgraph - reduced graph originating from a ring centers node
    graph - total graph of the Moltemplate file
    centers - a list of center nodes to be bonded
    start - node to find bonding ids for
    """
    bond_canidates = centers.copy()
    bond_canidates.remove(start)

    marked = []
    stack = [start]
    while stack:
        v = stack.pop()
        if v not in marked:
            # see if v is next to a ring center
            for w in rgraph[v]:
                if w in bond_canidates:
                    # remove ring from graph, return center id
                    marked += rgraph[v]
                    yield w
            # if not a ring, mark it and continue on
            marked.append(v)
            for w in rgraph[v]:
                if w not in marked:
                    stack.append(w)


def sort_bond_ids(id1, id2):
    if id1 < id2:
        first_atom = id1
        second_atom = id2
    else:
        first_atom = id2
        second_atom = id1
    return first_atom, second_atom


def add_center_atom(mol, ring_sizes=(5, 6, 7), method='centroid'):
    """
    Add a particle at the center of all rings defined by ring_sizes
    :param mol: molspace object
    :param ring_sizes: number of members in searched rings
    :param method: center particle location ('centroid' or 'CoM' - center of mass)
    :return:
    """
    centers = rings.get_rings(mol, ring_sizes=ring_sizes)

    # get atom parameters to copy, set indexes to to not overwrite
    atom_factory = mol.atoms.styles.atom_factory
    atom_id_iterator = max(mol.atoms) + 1
    u_type = max(mol.ff.masses) + 1
    mol.ff.masses[u_type] = FF.Parameters([20.18])  # neon dummy

    center_ids = []
    for ring in centers:
        # add a dummy atom denoting the center of a ring
        atom = atom_factory()
        # find coordinates of center
        if method == 'centroid':
            atom.x, atom.y, atom.z = ring.get_centroid()
        elif method == 'CoM':
            atom.x, atom.y, atom.z = ring.get_CoM()

        atom.type = u_type
        atom.id = atom_id_iterator
        center_ids.append((atom.id))  # add to list of ids for graph theory
        mol.atoms[atom_id_iterator] = atom  # add the atom to the system

        # add a bond from the dummy atom to each real atom in the ring
        for ring_atom in ring.members:
            id1, id2 = sort_bond_ids(ring_atom.id, atom_id_iterator)

            bond = mol.bonds.bond_factory()
            bond.ordered = [id1, id2]
            bond.type = 1
            mol.bonds.id_set(bond, id1, id2)

        atom_id_iterator = atom_id_iterator + 1

    return mol, center_ids


def coarsen_by_rings(mol, ring_sizes=(5, 6, 7), depth=1, center_method='centroid', remove_old_atoms=True):
    """
    :param mol: molspace object
    :param depth: depth of grid search to bond rings (1 is fused rings)
    :param center_method: method for centering rings ('centroid' or 'CoM')
    :param remove_old_atoms: remove starting atoms from molspace object, leaving only
        coarsened particles
    :return: molspace object with ring centers as particles
    """
    # 1 - Identify centers
    # 2 - Bond centers to their ring
    # 3 - Ensure they are a unique atom type
    mol, center_ids = add_center_atom(mol, ring_sizes=ring_sizes, method=center_method)

    # 4 - depth first search of atoms starting from centers
    graph = mol.generate_graph()
    for node in center_ids:
        rgraph = _reduced_graph(graph=graph, node=node, max_depth=depth)
        for chain_link in _dfs_chain(rgraph=rgraph, graph=graph, centers=center_ids, start=node):
            id1, id2 = sort_bond_ids(node, chain_link)

            # 5 - create a bond between dummy atoms
            # create a bond between centers
            bond = mol.bonds.bond_factory()
            bond.ordered = [id1, id2]
            bond.type = 1
            mol.bonds.id_set(bond, id1, id2)

        # 6 - clean up non-center atoms
        if remove_old_atoms:
            # remove non-center atoms:
            remove_ids = []
            for atom in mol.atoms:
                remove_ids.append(atom)

            for center in center_ids:
                remove_ids.remove(center)

            mol.remove_atoms(remove_ids)
            mol.update_elements()
    return mol

# -*- coding: utf-8 -*-
from collections import defaultdict
# from mooonpy.molspace.molspace import Molspace # cannot circular import. only need for types

class Cluster():
    # could inherit member as set?
    __slots__ = ['members', 'mass', 'composition']

    def __init__(self, in_set=None):
        if in_set is None:
            self.members = set()
        else:
            self.members = set(in_set)
        self.mass = 0.0
        self.composition = defaultdict(int)  # default 0 for default dict

    def __repr__(self):
        return 'Cluster object with {} members'.format(len(self.members))


class Clusters(dict):
    def __init__(self, molspace):
        self._molspace = molspace  # reference to parent? only use internally
        # dont compute here, molspace is mostly empty when initialized

    def __repr__(self):
        """Modify variable explorer label"""
        return 'Clusters class of {} Cluster objects'.format(len(self))

    def dict(self):
        """
        Return a dict of the cluster member sets
        """
        out_dict = {}
        for molid, cluster in self:
            out_dict[molid] = cluster.members
        return out_dict

    def mass_dict(self):
        """
        Returns dict of {molid:mass}
        """
        return {molid: cluster.mass for molid, cluster in self.items()}

    def assign_clusters(self):
        """
        Read molid from atom objects to populate clusters
        """
        self.clear()
        for _id, atom in self._molspace.atoms.items():
            molid = atom.molid
            if molid not in self:
                self[molid] = Cluster()
            self[molid].members.add(_id)
        return self

    def compute_clusters(self):
        """
        Compute molid from bond connectivity
        Note, this clears the object
        """
        neighbors = first_neighbor(self._molspace)
        loop_molid = 0
        lenatoms = len(self._molspace.atoms.keys())
        unsorted = defaultdict(set)
        checked = set()

        for id_, atom in self._molspace.atoms.items():
            if id_ not in checked:
                loop_molid += 1  # start at 1
                checked.add(id_)
                walking = [id_]  # depth 0 search
                unsorted[loop_molid].add(id_) # save start atom
                while len(walking) > 0:
                    to_walk = []
                    for w_id in walking:
                        for n_id in neighbors[w_id]:
                            if n_id not in checked:
                                checked.add(n_id)
                                unsorted[loop_molid].add(n_id)
                                to_walk.append(n_id)  # N+1 depth search
                    walking = to_walk  # next depth of search
                if len(checked) == lenatoms:
                    break  # exit if all clusters explored

        # sorted_molid = {molid: unsorted[molid] for molid in
        #                 sorted(unsorted, key=lambda molid: len(unsorted[molid]), reverse=True)}

        sorted_molid = dict(sorted(unsorted.items(), key=lambda item: len(item[1]),reverse=True))

        # Note this is not deterministic if multiple clusters of the same size exist

        self.clear()
        new_molid = 0
        for molid, members in sorted_molid.items():
            new_molid +=1
            self[new_molid] = Cluster(members)
            ## update per-atom molid
            for id_ in members:
                self._molspace.atoms[id_].molid = new_molid

        return self

    def compute_info(self):
        """
        Compute the count of atom types in each cluster, and total mass
        """
        try:
            mass_lookup = {type_: param.coeffs[0] for type_, param in self._molspace.ff.masses.items()}
        except:
            mass_lookup = None
        for molid, cluster in self.items():
            cluster.composition.clear()
            for id_ in cluster.members:
                type_ = self._molspace.atoms[id_].type
                cluster.composition[type_] += 1

            cluster.mass = 0.0
            for type_, count in cluster.composition.items():
                if mass_lookup:  # skips if not assigned
                    cluster.mass += mass_lookup[type_] * count


def first_neighbor(molspace):
    """
    Function to compute single step per-atom graph
    """
    # neighbors = defaultdict(list) # or set? gor for ease not speed here
    neighbors = {}
    for id_, atom in molspace.atoms.items():  # need loop here for lone atoms
        # neighbors[id_] = []
        neighbors[id_] = set()

    bond = None  # dummy
    try:
        for bond in molspace.bonds.keys():
            # neighbors[bond[1]].append(bond[0])
            # neighbors[bond[0]].append(bond[1])
            neighbors[bond[1]].add(bond[0])
            neighbors[bond[0]].add(bond[1])
    except:
        raise Exception(f'Bond {bond} does not have associated atoms')

    return neighbors

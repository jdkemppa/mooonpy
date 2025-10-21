import mooonpy.molspace.molspace
import numpy as np

class Ring(object):
    """
    Object containing molfile info for a simulation box as well as
    a dictionary contining atom objects for each ring
    """

    def __init__(self, in_set=None, mol=None):
        self.mol = mol
        if in_set is None:
            self.members = set()
        else:
            self.members = set(in_set)

    def get_CoM(self):
        self.reset_ring_image()

        num = (0, 0, 0)
        ring_mass = 0
        for atom in self.members:
            mass = mol.ff.masses[atom.type].coeffs[0]
            atom_coords = (atom.x, atom.y, atom.z)
            num = tuple(num[i] + atom_coords[i] * mass for i in range(3))

            ring_mass = ring_mass + mass
        return tuple(num[i] / ring_mass for i in range(3))

    def get_centroid(self):
        self.reset_ring_image()
        num = (0, 0, 0)
        for atom in self.members:
            atom_coords = (atom.x, atom.y, atom.z)
            num = tuple(num[i] + atom_coords[i] for i in range(3))

        return tuple(num[i] / len(self.members) for i in range(3))

    def reset_ring_image(self):
        positions = np.empty((len(self.members), 3,))
        for i, atom in enumerate(self.members):
            positions[i, :] = [atom.x, atom.y, atom.z]
        ref = positions[0, :].copy()
        positions -= ref  # compute w.r.t lowest atom ID

        box = self.mol.atoms.box
        lx = box.xhi - box.xlo
        ly = box.yhi - box.ylo
        lz = box.zhi - box.zlo

        for coord in positions:
            for i, (c, l) in enumerate(zip(coord, [lx, ly, lz])):
                if c > l / 2:
                    coord[i] -= l
                elif c < -l / 2:
                    coord[i] += l

        positions += ref  # un-zero atom[0] coords and adjust other atoms accordingly

        # apply corrected coordinates to atoms in ring
        for atom, newpos in zip(self.members, positions):
            atom.x = newpos[0]
            atom.y = newpos[1]
            atom.z = newpos[2]


def get_rings(mol, ring_sizes=(5, 6, 7)):
    """
    Returns a list of Ring objects
    :param mol: Molspace file with bonds
    :param ring_sizes: N-membered rings to find
    :return: list of Ring objects
    """
    rings_list = mol.find_rings(ring_sizes=ring_sizes)
    rings = []
    for ringIDs in rings_list:
        # ring is a list of atom objects
        ring = Ring({mol.atoms[a] for a in list(ringIDs)}, mol=mol)
        ring.reset_ring_image()
        rings.append(ring)
    return rings

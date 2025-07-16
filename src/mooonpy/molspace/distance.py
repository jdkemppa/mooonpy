# -*- coding: utf-8 -*-
from collections import defaultdict
from dataclasses import dataclass

from mooonpy.molspace.topology import Bonds, Angles, Dihedrals, Impropers
from mooonpy.molspace.atoms import Atoms

# from ..tools.math_utils import MixingRule
import math
import numpy as np


@dataclass
class Domain(list):
    """
    Class to contain atoms within a domain as a list of atom IDs
    """

    def __init__(self, atom_list=None):
        if atom_list is None:
            atom_list = []
        super(Domain, self).__init__(atom_list)
        self.x: float = 0.0
        self.y: float = 0.0
        self.z: float = 0.0
        self.neighbors: dict = {}


@dataclass
class Pair():
    def __init__(self, dx, dy, dz, distance):
        self.dx: float = dx
        self.dy: float = dy
        self.dz: float = dz
        self.distance: float = distance


class Pairs(dict):
    def __init__(self, from_dict=None, to_dict=None):
        if from_dict is None:
            super(Pairs, self).__init__()
        else:
            super(Pairs, self).__init__(from_dict)

    def update_bonds(self, bonds, vect=True, ignore_missing=False):
        """
        Update Bonds with length attribute computed from pair interactions

        :param bonds: Bonds object to update
        :type bonds: Bonds
        :param vect: If True, updates bond.vect, may be disabled for speed
        :type vect: bool
        """
        for key, bond in bonds.items():
            try:
                pair = self[key]
                bond.dist = pair.distance
                if vect:
                    bond.vect = (pair.dx, pair.dy, pair.dz)
            except:
                if not ignore_missing:
                    raise KeyError('Bond has no matching key in Pairs, length exceeded or it may not exist')

    def filter_cutoff(self, atoms=None, bonds=None, cutoff=None, mode=None):
        """
        Return modified Pairs list after rule, may also update atoms or bonds
        """

        if not isinstance(atoms, Atoms):
            raise TypeError('atoms must be a Atoms object')


def pairs_from_bonds(atoms, bonds, periodicity='ppp'):
    """
    Compute pairs from bonds using minimum image convention
    If bonds span more than half the box span in a periodic direction, the bond vector
    """
    if not isinstance(atoms, Atoms):
        raise TypeError('atoms must be a Atoms object')
    if not isinstance(bonds, Bonds):
        raise TypeError('bonds must be a Bond object')

    ## There may be a quicker way to do this, or one that works periodically for small cells
    box = atoms.box
    h, h_inv, boxlo, boxhi = box.get_transformation_matrix()
    lx, ly, lz = box.get_lengths()

    fractionals = {}
    for id_, atom in atoms.items():
        fractionals[id_] = box.pos2frac(atom.x, atom.y, atom.z, h_inv, boxlo)

    pairs = Pairs()
    for key, bond in bonds.items():
        try:
            f_A = fractionals[key[0]]
            f_B = fractionals[key[1]]
        except:
            raise KeyError(f'Bond key {key} has no matching key in Atoms')
        du_x = f_B[0] - f_A[0]
        du_y = f_B[1] - f_A[1]
        du_z = f_B[2] - f_A[2]
        if periodicity[0] == 'p':
            if du_x > 0.5:
                du_x -= 1
            elif du_x < -0.5:
                du_x += 1
        if periodicity[1] == 'p':
            if du_y > 0.5:
                du_y -= 1
            elif du_y < -0.5:
                du_y += 1
        if periodicity[2] == 'p':
            if du_z > 0.5:
                du_z -= 1
            elif du_z < -0.5:
                du_z += 1

        ## Transform fractional vector. not using function because no boxlo
        dx = h[0] * du_x + h[5] * du_y + h[4] * du_z
        dy = h[1] * du_y + h[3] * du_z
        dz = h[2] * du_z

        distance2 = dx * dx + dy * dy + dz * dz
        pairs[key] = Pair(dx, dy, dz, math.sqrt(distance2))

    pairs.update_bonds(bonds)  # update bond object
    return pairs


def pairs_from_domains(atoms, cutoff, domains, fractionals):
    if not isinstance(atoms, Atoms):
        raise TypeError('atoms must be a Atoms object')
    box = atoms.box
    h, h_inv, boxlo, boxhi = box.get_transformation_matrix()

    # pairs = {}
    pairs = Pairs()
    cutoff_pow2 = cutoff * cutoff

    for box_index, domain in domains.items():  # this can be parallelized
        for neighbor_index, image_shift in domain.neighbors.items():
            for atom_A in domain:  # in list
                ## setup 1 A atom at a time
                if image_shift == (0, 0, 0):  # no shift, get positions skipping box steps
                    Atom_A = atoms[atom_A]  # get class instance
                    pos_x = Atom_A.x
                    pos_y = Atom_A.y
                    pos_z = Atom_A.z
                else:  # get positions after shift
                    ux, uy, uz = fractionals[atom_A]
                    ux += image_shift[0]
                    uy += image_shift[1]
                    uz += image_shift[2]
                    pos_x, pos_y, pos_z = box.frac2pos(ux, uy, uz, h, boxlo)
                for atom_B in domains[neighbor_index]:
                    Atom_B = atoms[atom_B]  # get class instance
                    if atom_A == atom_B:
                        continue
                    elif atom_A < atom_B:  # from low to high
                        dx = Atom_B.x - pos_x
                        dy = Atom_B.y - pos_y
                        dz = Atom_B.z - pos_z
                        key = (atom_A, atom_B)
                    else:  # reverse vector
                        dx = pos_x - Atom_B.x
                        dy = pos_y - Atom_B.y
                        dz = pos_z - Atom_B.z
                        key = (atom_B, atom_A)
                    distance2 = dx * dx + dy * dy + dz * dz
                    if distance2 < cutoff_pow2:
                        pairs[key] = Pair(dx, dy, dz, math.sqrt(distance2))
    return pairs


def domain_decomp_13(atoms, cutoff, whitelist=None, blacklist=None, periodicity='ppp'):
    """
    Setup domains for pairwise distance computation. This uses a 3x3x3 grid, where each domain checks self
    and half of the 26 adjacent domains, hence 13 others. Checking order priority is -z, -x, -y (make image eventually)
    This algorithm is optimized for short cutoffs under ~5 angstroms. Longer cutoffs would be faster with a
    DD_62 algorithm similar to https://docs.lammps.org/Developer_par_neigh.html


    ..warning :: Atoms must be correctly wrapped before calling this function.
        Highly skewed triclinic does not currently set the minimum domain sizes correctly

    ..TODO :: Validate triclinic minimum domains

    """
    if not isinstance(atoms, Atoms):
        raise TypeError('atoms must be a Atoms object')
    if not isinstance(periodicity, str) or len(periodicity) != 3:
        raise TypeError('periodicity must be a string of form "pp?"')

    # atoms = atoms.copy()
    # atoms.wrap() # user should use these externally first

    ## Setup fractional conversions and domain sizes
    box = atoms.box
    h, h_inv, boxlo, boxhi = box.get_transformation_matrix()
    lx, ly, lz = box.get_lengths()

    nx = int(lx // cutoff)  # number of domains
    ny = int(ly // cutoff)
    nz = int(lz // cutoff)
    if nx == 0: nx = 1
    if ny == 0: ny = 1
    if nz == 0: nz = 1

    ## number of domains for skewed triclinic, needs testing and 3D checks
    # nx = int(lx // (cutoff / math.cos(math.atan(box.xy / ly)) / math.cos(math.atan(box.xz / lz))))
    # ny = int(ly // (cutoff / math.cos(math.atan(box.yz / lz))))
    # nz = int(lz // cutoff)

    if nx < 3 and periodicity[0] == 'p':
        raise Exception(f'cutoff {cutoff} is too large for box x {lx} to create 3 periodic domains in the x direction')
    if ny < 3 and periodicity[1] == 'p':
        raise Exception(f'cutoff {cutoff} is too large for box y {ly} to create 3 periodic domains in the y direction')
    if nz < 3 and periodicity[2] == 'p':
        raise Exception(f'cutoff {cutoff} is too large for box z {lz} to create 3 periodic domains in the z direction')

    fx = 1 / nx  # use floor div later to get domain index
    fy = 1 / ny
    fz = 1 / nz

    ## fractionalize atoms and assign to domains
    fractionals = {}
    domains = defaultdict(Domain)

    for id_, atom in atoms.items():
        if whitelist is not None and id_ not in whitelist:
            continue
        elif blacklist is not None and id_ in blacklist:
            continue

        frac = box.pos2frac(atom.x, atom.y, atom.z, h_inv, boxlo)
        fractionals[id_] = frac
        dx = int(frac[0] // fx)  # box index [0,nx) exclusive if wrapped correctly
        dy = int(frac[1] // fy)
        dz = int(frac[2] // fz)
        box_index = (dx, dy, dz)
        if min(box_index) < 0 or dx >= nx or dy >= ny or dz >= nz:
            raise Exception('Domain index is out of bounds, use Atoms.wrap() before computation')

        domains[box_index].append(id_)  # this is the only problem spot if parallelized

    # self+13 neighboring domains. down, east, south priority ordering.
    neighbor_shifts = [(0, 0, 0), (0, 0, -1), (0, -1, 0), (0, -1, -1), (-1, 0, 0), (-1, 0, -1), (-1, -1, 0),
                       (-1, -1, -1), (-1, 1, 0), (-1, 1, -1), (0, 1, -1), (1, 1, -1), (1, 0, -1), (1, -1, -1)]

    for box_index, domain in domains.items():
        # not sure about use case for these, but these are the lower corner positions
        domain.x, domain.y, domain.z = box.frac2pos(box_index[0] / nx, box_index[1] / ny, box_index[2] / nz, h, boxlo)

        ## Setup neighbor domains to loop through
        for neighbor_shift in neighbor_shifts:
            neighbor_index = []
            image_shift = []
            ## loop directions
            for self_i, shift_i, n_i, period_i in zip(box_index, neighbor_shift, (nx, ny, nz), periodicity):
                shifted = self_i + shift_i
                if shifted >= 0 and shifted < n_i:  # middle of box or looking in
                    neighbor_index.append(shifted)
                    image_shift.append(0)

                elif period_i != 'p':
                    break  # non-periodic, so this neighbor is skipped

                elif shifted < 0:
                    neighbor_index.append(n_i - 1)  # go to top
                    image_shift.append(1)  # add 1 box length in this direction when comparing

                else:  # shifted >= n_i:
                    neighbor_index.append(0)  # go to bottom
                    image_shift.append(-1)  # subtract 1 box length
                # else:
                #     raise Exception(f'Something went wrong, this should be unreachable')

            else:  # Skipped by break if non-periodic
                neighbor_index = tuple(neighbor_index)
                if neighbor_index in domains:  # skip empty domains later, faster for sparse boxes
                    # save to this domain neighbors
                    domain.neighbors[neighbor_index] = tuple(image_shift)

    return domains, fractionals


def ADI_from_bonds(bonds, angles=None, dihedrals=None, impropers=None):
    if not isinstance(bonds, Bonds):
        raise TypeError('bonds must be a Bond object')

    vectors = {}
    for key, bond in bonds.items():
        if not hasattr(bond, 'vect'):
            raise Exception('bond has no computed vector, use "compute_bond_length" first')
        vectors[key] = np.array(bond.vect) / bond.dist

    if isinstance(angles, Angles):
        for key, angle in angles.items():
            b1 = bonds.key_rule(angle.ordered[1], angle.ordered[0])  # outward if keys match
            b2 = bonds.key_rule(angle.ordered[1], angle.ordered[2])  # may break if rules change. This will need tests

            v1 = vectors[b1]
            v2 = vectors[b2]
            if b1[0] != angle.ordered[1]: v1 = -v1  # if tip is center mirror a copy
            if b2[0] != angle.ordered[1]: v2 = -v2  # *= changes in the dict

            dot = np.dot(v1, v2) # / bonds[b1].dist / bonds[b2].dist
            angle.theta = np.arccos(dot) * 180 / np.pi
            ## could do a cross product to get normal?
    else:
        pass

    if isinstance(dihedrals, Dihedrals):
        for key, dihedral in dihedrals.items():
            b1 = bonds.key_rule(dihedral.ordered[0], dihedral.ordered[1])
            b2 = bonds.key_rule(dihedral.ordered[1], dihedral.ordered[2])
            b3 = bonds.key_rule(dihedral.ordered[2], dihedral.ordered[3])

            v1 = vectors[b1]
            v2 = vectors[b2]
            v3 = vectors[b3]
            if b1[0] != dihedral.ordered[1]: v1 = -v1  # out from 1
            if b2[0] != dihedral.ordered[1]: v2 = -v2  # out from 1
            if b3[0] != dihedral.ordered[2]: v3 = -v3  # out from 2

            vl = np.cross(v2, v1)  # I think these give the correct sign
            vr = np.cross(v2, v3)
            matrix = np.vstack([vl, vr, v2])
            triple_product = np.linalg.det(matrix)
            # triple_product = (vl[0] * (vr[1] * v2[2] - vr[2] * v2[1]) -
            #                   vl[1] * (vr[0] * v2[2] - vr[2] * v2[0]) +
            #                   vl[2] * (vr[0] * v2[1] - vr[1] * v2[0])) # slower

            cos_theta = np.dot(vl, vr)  # not normalized, cancels out
            sin_theta = triple_product # / bonds[b2].dist # inv vl and inv vr would also be here
            a = np.atan2(sin_theta, cos_theta) * 180 / np.pi

            dihedral.phi = a
    else:
        pass

    if isinstance(impropers, Impropers):
        for key, improper in impropers.items():
            b1 = bonds.key_rule(improper.ordered[1], improper.ordered[0])
            b2 = bonds.key_rule(improper.ordered[1], improper.ordered[2])
            b3 = bonds.key_rule(improper.ordered[1], improper.ordered[3])

            v1 = vectors[b1]
            v2 = vectors[b2]
            v3 = vectors[b3]

            if b1[0] != improper.ordered[1]: v1 = -v1  # all center out
            if b2[0] != improper.ordered[1]: v2 = -v2
            if b3[0] != improper.ordered[1]: v3 = -v3

            ## Slow
            norm1 = np.cross(v2, v3)
            norm1 /= np.linalg.norm(norm1)
            norm2 = np.cross(v3, v1) # order matches LAMMPS class2.cpp
            norm2 /= np.linalg.norm(norm2)
            norm3 = np.cross(v1, v2)
            norm3 /= np.linalg.norm(norm3)

            ch1 = np.arcsin(np.dot(norm1, v1)) # / bonds[b1].dist)
            ch2 = np.arcsin(np.dot(norm2, v2)) # / bonds[b2].dist)
            ch3 = np.arcsin(np.dot(norm3, v3)) # / bonds[b3].dist)

            # ch1 = np.arcsin(np.linalg.det(np.vstack([v2, v3, v1])) / (bonds[b1].dist*bonds[b2].dist*bonds[b3].dist))
            # ch2 = np.arcsin(np.linalg.det(np.vstack([v3, v1, v2])) / (bonds[b1].dist*bonds[b2].dist*bonds[b3].dist))
            # ch3 = np.arcsin(np.linalg.det(np.vstack([v1, v2, v3])) / (bonds[b1].dist*bonds[b2].dist*bonds[b3].dist))
            ## not work because norm1 mag is not the same as norm of v2 and v3



            improper.chi = (ch1 + ch2 + ch3) * 180 / (3 * np.pi)
    else:
        pass


def BADI_by_type(mol, type_label=False, comp_bond=True, comp_angle=True, comp_dihedral=True, comp_improper=True):
    if comp_bond:
        bond_types = mol.ff.bond_coeffs
        if type_label:
            bond_hist = {coeff.type_label: [] for coeff in bond_types.values()}
        else:
            bond_hist = {type_: [] for type_ in bond_types.keys()}

        for bond in mol.bonds.values():
            key = bond.type
            if type_label: key = bond_types[key].type_label
            bond_hist[key].append(bond.dist)
    else:
        bond_hist = {}

    if comp_angle:
        angle_types = mol.ff.angle_coeffs
        if type_label:
            angle_hist = {coeff.type_label: [] for coeff in angle_types.values()}
        else:
            angle_hist = {type_: [] for type_ in angle_types.keys()}

        for angle in mol.angles.values():
            key = angle.type
            if type_label: key = angle_types[key].type_label
            angle_hist[key].append(angle.theta)
    else:
        angle_hist = {}

    if comp_dihedral:
        dihedral_types = mol.ff.dihedral_coeffs
        if type_label:
            dihedral_hist = dihedral_hist = {coeff.type_label: [] for coeff in dihedral_types.values()}
        else:
            dihedral_hist = {type_: [] for type_ in dihedral_types.keys()}

        for dihedral in mol.dihedrals.values():
            key = dihedral.type
            if type_label: key = dihedral_types[key].type_label
            dihedral_hist[key].append(dihedral.phi)
    else:
        dihedral_hist = {}

    if comp_improper:
        improper_types = mol.ff.improper_coeffs
        if type_label:
            improper_hist = {coeff.type_label: [] for coeff in improper_types.values()}
        else:
            improper_hist = {type_: [] for type_ in improper_types.keys()}

        for improper in mol.impropers.values():
            key = improper.type
            if type_label: key = improper_types[key].type_label
            improper_hist[key].append(improper.chi)
    else:
        improper_hist = {}

    return bond_hist, angle_hist, dihedral_hist, improper_hist

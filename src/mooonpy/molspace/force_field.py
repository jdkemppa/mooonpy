# -*- coding: utf-8 -*-
import mooonpy.molspace.periodic_table
import numpy as np
from copy import deepcopy


class Parameters(object):
    def __init__(self, coeffs):
        self.coeffs: list = coeffs
        self.style: str = ''
        self.comment: str = ''
        self.type_label: str = ''
        self.element: str = ''


class Coefficients(dict):
    def __init__(self, keyword, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.style: str = ''
        self.keyword = keyword


class ForceField(object):
    def __init__(self, **kwargs):
        # Build this object with some composition
        self.masses: Coefficients = Coefficients('Masses')  # {type : Parameters-object}
        self.pair_coeffs: Coefficients = Coefficients('Pair Coeffs')  # {type : Parameters-object}

        self.bond_coeffs: Coefficients = Coefficients('Bond Coeffs')  # {type : Parameters-object}
        self.angle_coeffs: Coefficients = Coefficients('Angle Coeffs')  # {type : Parameters-object}
        self.dihedral_coeffs: Coefficients = Coefficients('Dihedral Coeffs')  # {type : Parameters-object}
        self.improper_coeffs: Coefficients = Coefficients('Improper Coeffs')  # {type : Parameters-object}

        self.bondbond_coeffs: Coefficients = Coefficients('BondBond Coeffs')  # {type : Parameters-object}
        self.bondangle_coeffs: Coefficients = Coefficients('BondAngle Coeffs')  # {type : Parameters-object}
        self.angleangletorsion_coeffs: Coefficients = Coefficients(
            'AngleAngleTorsion Coeffs')  # {type : Parameters-object}
        self.endbondtorsion_coeffs: Coefficients = Coefficients('EndBondTorsion Coeffs')  # {type : Parameters-object}
        self.middlebondtorsion_coeffs: Coefficients = Coefficients(
            'MiddleBondTorsion Coeffs')  # {type : Parameters-object}
        self.bondbond13_coeffs: Coefficients = Coefficients('BondBond13 Coeffs')  # {type : Parameters-object}
        self.angletorsion_coeffs: Coefficients = Coefficients('AngleTorsion Coeffs')  # {type : Parameters-object}
        self.angleangle_coeffs: Coefficients = Coefficients('AngleAngle Coeffs')  # {type : Parameters-object}

        # Boolean to check if type labels have been read-in
        self.has_type_labels = False

    def coeffs_factory(self, coeffs=None):
        if coeffs is None: coeffs = []
        return Parameters(coeffs)

    def get_per_line_styles(self, coeff):
        lines = {}  # {'TypeID':'per-line-style'}
        potential = getattr(self, coeff)
        for i in potential:
            parameters = potential[i]
            lines[i] = parameters.style
        return lines

    def copy(self):
        return deepcopy(self)

    def elements2mass(self, atoms, elements=None, change_types=False):
        """
        Update masses Coefficients using per-atom element labels
        :param Atoms: Input Atoms object to loop through
        :param elements: list of elements to be used (1 indexed in dict output), if None, will loop through atoms and use periodic_table.py ordering
        :param change_types: Bool, will update each Atom.type with the new type for each element
        """
        from mooonpy.molspace.atoms import Atoms
        if not isinstance(atoms, Atoms):
            raise TypeError('atoms must be a Atoms object')
        pt = mooonpy.molspace.periodic_table.Elements()
        if elements is None:
            elements_used = set()
            for atom in atoms.values():
                elements_used.add(atom.element)
            elements = []
            for ele, Element in pt.elements.items():
                if ele in elements_used:
                    elements.append(ele)
        # print(elements)
        masses = self.masses
        masses.clear()
        ii = 0
        type_map = {}
        for ele in elements:
            ii += 1
            coeff = self.coeffs_factory()
            coeff.element = ele
            coeff.coeffs = [pt.elements[ele].masses[0]]
            masses[ii] = coeff
            type_map[ele] = ii
        if change_types:
            for id_, atom in atoms.items():
                atom.type = type_map[atom.element]

    def reax(self, elements=None, dummybond=False):
        if elements is None:
            elements = {}
            for type_, param in self.masses.items():
                ele = param.element
                if ele not in elements:
                    elements[ele] = param.coeffs

        out_ff = ForceField()
        masses = out_ff.masses
        for ii, element in enumerate(elements):
            ii += 1
            masses[ii] = Parameters(elements[element])
            masses[ii].element = ele

        if dummybond:
            out_ff.bond_coeffs[1] = self.coeffs_factory()

        return out_ff


## these may move somewhere else? or be callable from coeff class directly?
def class2_dihedral(coeffs, angles):
    energy = np.zeros_like(angles)
    for ii in range(3):
        if coeffs[ii * 2] != 0:
            energy += coeffs[ii * 2] * (1 - np.cos(np.deg2rad((ii + 1) * angles - coeffs[ii * 2 + 1])))
    return energy

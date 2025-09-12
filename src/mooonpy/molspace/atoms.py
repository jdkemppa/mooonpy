# -*- coding: utf-8 -*-
"""
This module provides a class to organize atoms information
"""
from copy import deepcopy

from .box import Box
from .atom_styles import Styles

from typing import Dict, Tuple


class Atoms(dict):
    def __init__(self, astyles, **kwargs):
        super().__init__(**kwargs)
        self.style: str = 'full'  # will default to full and update when needed

        # Build this object with some composition
        self.box: Box = Box()
        self.styles: Styles = Styles(astyles)

    def shift(self, sx=0, sy=0, sz=0):
        return

    def copy(self):
        return deepcopy(self)

    def wrap(self):
        """
        Wrap all coordinates so all atoms are inside the box, and indexes box image appropriately.

        .. seealso:: box.get_transformation_matrix, box.pos2frac, box.frac2pos

        ..TODO::
            This should work for triclinic, but is currently untested
        """
        h, h_inv, boxlo, boxhi = self.box.get_transformation_matrix()
        for id_, atom in self.items():
            ux, uy, uz = self.box.pos2frac(atom.x, atom.y, atom.z, h_inv, boxlo)
            atom.ix += int(ux // 1)
            atom.iy += int(uy // 1)
            atom.iz += int(uz // 1)

            pos_x, pos_y, pos_z = self.box.frac2pos(ux % 1, uy % 1, uz % 1, h, boxlo)
            atom.x = pos_x
            atom.y = pos_y
            atom.z = pos_z
        return

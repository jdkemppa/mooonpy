# -*- coding: utf-8 -*-
"""
This module provides a class to organize atoms information
"""
import warnings
from copy import deepcopy

from .box import Box
from .atom_styles import Styles

from typing import Dict, Iterable, Optional, Tuple


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

    def _select(self, atom_ids):
        """
        Resolve a per-call atom-id selection. ``'all'`` (or ``None``) means every
        atom. Otherwise iterate the supplied ids.
        """
        if atom_ids is None or atom_ids == 'all':
            return self.values()
        return (self[id_] for id_ in atom_ids)

    def centroid(self, atom_ids='all') -> Tuple[float, float, float]:
        """
        Unweighted centroid of the selected atoms.

        :param atom_ids: ``'all'`` or an iterable of atom ids.
        :return: ``(x, y, z)`` centroid.
        """
        n = 0
        cx = cy = cz = 0.0
        for atom in self._select(atom_ids):
            cx += atom.x
            cy += atom.y
            cz += atom.z
            n += 1
        if n == 0:
            return 0.0, 0.0, 0.0
        return cx / n, cy / n, cz / n

    def COM(self, masses, atom_ids='all') -> Tuple[float, float, float]:
        """
        Mass-weighted center of mass of the selected atoms.

        :param masses: A ``{type: mass}`` dict. Pass either a plain dict or
                       ``mol.ff.masses`` (the ``Parameters`` object's
                       ``coeffs[0]`` is unwrapped automatically).
        :param atom_ids: ``'all'`` or an iterable of atom ids.
        :return: ``(x, y, z)`` center of mass.
        :raises KeyError: If a selected atom's ``type`` is not in ``masses``.
        """
        # Normalize masses input: accept ff.masses-style {type: Parameters} or plain {type: float}.
        mass_lookup = {}
        for type_, val in masses.items():
            if hasattr(val, 'coeffs'):
                mass_lookup[type_] = val.coeffs[0]
            else:
                mass_lookup[type_] = val

        total = 0.0
        mx = my = mz = 0.0
        for atom in self._select(atom_ids):
            m = mass_lookup[atom.type]
            total += m
            mx += m * atom.x
            my += m * atom.y
            mz += m * atom.z
        if total == 0.0:
            return 0.0, 0.0, 0.0
        return mx / total, my / total, mz / total

    def move(self, vect=None, mode='offset', atom_ids='all'):
        """
        Translate or scale atom coordinates in place.

        :param vect: Displacement / target / scale vector. ``None`` is a no-op.
        :type vect: tuple[float, float, float] or None
        :param mode: How ``vect`` is interpreted:

            * ``'offset'``   — add ``vect`` to each atom's coordinates.
            * ``'centroid'`` — shift the selection so its centroid lands at
              ``vect`` (computed via :meth:`centroid`).
            * ``'scale'``    — multiply each coordinate component by ``vect``.
              This also scales the simulation :class:`Box` parameters about the
              origin (``xlo *= vx`` etc., tilt factors scale with the
              relevant axis).

        :param atom_ids: ``'all'`` or an iterable of atom ids to apply the move
                         to. ``'scale'`` always scales the box regardless of
                         the selection.

        For mass-weighted re-centering, compute the COM externally with
        :meth:`COM` and pass the desired delta with ``mode='offset'``.
        """
        if vect is None:
            return
        vx, vy, vz = float(vect[0]), float(vect[1]), float(vect[2])

        if mode == 'offset':
            for atom in self._select(atom_ids):
                atom.x += vx
                atom.y += vy
                atom.z += vz
            return

        if mode == 'centroid':
            cx, cy, cz = self.centroid(atom_ids)
            dx, dy, dz = vx - cx, vy - cy, vz - cz
            for atom in self._select(atom_ids):
                atom.x += dx
                atom.y += dy
                atom.z += dz
            return

        if mode == 'scale':
            for atom in self._select(atom_ids):
                atom.x *= vx
                atom.y *= vy
                atom.z *= vz
            # Box scales about the origin. Tilt factors scale with the axis
            # they parameterize: xy and xz are x-direction offsets (scale by
            # vx), yz is a y-direction offset (scale by vy).
            self.box.xlo *= vx
            self.box.xhi *= vx
            self.box.ylo *= vy
            self.box.yhi *= vy
            self.box.zlo *= vz
            self.box.zhi *= vz
            self.box.xy *= vx
            self.box.xz *= vx
            self.box.yz *= vy
            return

        raise ValueError(
            f"Unknown move mode {mode!r}; expected 'offset', 'centroid', or 'scale'"
        )

    def shrinkwrap(self, pad=0.0):
        """
        Resize :attr:`box` to tightly enclose every atom, with optional padding
        on each face.

        :param pad: Distance (Å) added to each face beyond the atom extent.
        :type pad: float
        :raises ValueError: If the box is triclinic (any of ``xy``/``xz``/``yz``
                            is nonzero); shrinkwrap on tilted boxes is
                            ill-defined here.
        """
        box = self.box
        if box.xy or box.xz or box.yz:
            raise ValueError(
                "shrinkwrap requires an orthorhombic box; got tilt factors "
                f"xy={box.xy}, xz={box.xz}, yz={box.yz}"
            )
        if not self:
            return
        xs = [a.x for a in self.values()]
        ys = [a.y for a in self.values()]
        zs = [a.z for a in self.values()]
        box.xlo = min(xs) - pad
        box.xhi = max(xs) + pad
        box.ylo = min(ys) - pad
        box.yhi = max(ys) + pad
        box.zlo = min(zs) - pad
        box.zhi = max(zs) + pad

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

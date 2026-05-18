# -*- coding: utf-8 -*-
"""
Reaction-template tools for mooonpy. Builds pre/post Molspaces and a map data
structure for use with LAMMPS ``fix bond/react`` (the REACTER workflow).
"""
import importlib

__all__ = ['fragment',
           'reaction',
           'mapfile',
           'ntafile',
]

for name in __all__:
    module = importlib.import_module(f'.{name}', __package__)
    globals()[name] = module

from mooonpy.template.fragment import Fragment
from mooonpy.template.reaction import ReactionTemplate
from mooonpy.template.mapfile import write_mapfile
from mooonpy.template.ntafile import write_ntafile, canonicalize_ntas

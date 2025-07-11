# -*- coding: utf-8 -*-
import importlib

__all__ = ['graph',
]

for name in __all__:
    module = importlib.import_module(f'.{name}', __package__)
    globals()[name] = module
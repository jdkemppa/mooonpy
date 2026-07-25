# -*- coding: utf-8 -*-
import importlib
from typing import TYPE_CHECKING

_aliases = {
    'Molspace': '.molspace.molspace',
    'DocExamples': '.molspace',
    'Thermospace': '.thermospace.thermospace',
    'Path': '.tools.file_utils',
    'ReactionTemplate': '.template.reaction',
    # submodules
    'molspace': '.molspace',
    'programs': '.programs',
    'thermospace': '.thermospace',
    'tools': '.tools',
    'xrdspace': '.xrdspace',
    'fitting': '.fitting',
    'template': '.template',
}

__all__ = list(_aliases.keys())

if TYPE_CHECKING: # style/type hints
    from .fitting import fitting as fitting
    from .molspace import doc_examples as DocExamples
    from .molspace.molspace import Molspace
    from .programs import programs as programs
    from .template.reaction import ReactionTemplate
    from .thermospace.thermospace import Thermospace
    from .tools.file_utils import Path


def __getattr__(name: str):
    if name in _aliases:
        module = importlib.import_module(_aliases[name], __package__)
        try:
            obj = getattr(module, name)
        except AttributeError:
            if module.__name__.split('.')[-1] == name:
                obj = module
            else:
                raise

        globals()[name] = obj # cache module for later lookups
        return obj

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return sorted(list(globals().keys()) + __all__)
# -*- coding: utf-8 -*-
import importlib
import os


from .rcsetup import rcParams

# Generate "aliased" imports
from .molspace.molspace import Molspace as Molspace
from .molspace import doc_examples as DocExamples

from .thermospace.thermospace import Thermospace as Thermospace ## TDM
from .tools.file_utils import Path as Path ## TDM


from .rcsetup import rcParams

__all__ = ['Molspace',
           'DocExamples',
           'Thermospace',
           'Path',
]


# Get the current directory (the package root)
_package_dir = os.path.dirname(__file__)


# List of submodules/folders to import
_submodules = ['guis', 'molspace', 'programs', 'thermospace', 'tools', 'xrdspace','fitting']
for name in _submodules:
    if os.path.isdir(os.path.join(_package_dir, name)):
        module = importlib.import_module(f'.{name}', __package__)
        globals()[name] = module
        __all__.append(name)
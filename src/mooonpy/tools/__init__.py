# -*- coding: utf-8 -*-
"""
Kept as lazy imports since some of the math and signals analysis can bog
down import times
"""

import importlib
from typing import TYPE_CHECKING

_utils = {
    'Path': '.file_utils',
    'Cartesian': '.loop_utils',
    'ColorMap': '.misc_utils',
    'Table': '.tables',
    # Submodules
    'file_utils': '.file_utils',
    'loop_utils': '.loop_utils',
    'math_utils': '.math_utils',
    'misc_utils': '.misc_utils',
    'signals': '.signals',
    'string_utils': '.string_utils',
    'tables': '.tables',
}

__all__ = list(_utils.keys())

# PEP 484: IDE Autocomplete and Type Checking
if TYPE_CHECKING:
    from .file_utils import Path
    from .loop_utils import Cartesian
    from .misc_utils import ColorMap
    from .tables import Table

    from . import file_utils as file_utils
    from . import loop_utils as loop_utils
    from . import math_utils as math_utils
    from . import misc_utils as misc_utils
    from . import signals as signals
    from . import string_utils as string_utils
    from . import tables as tables

# PEP 562: Dynamic Attribute Access
def __getattr__(name: str):
    if name in _utils:
        module = importlib.import_module(_utils[name], __package__)
        try:
            # Try to grab the specific class/function/variable from the module
            obj = getattr(module, name)
        except AttributeError:
            # Fallback: Handles importing submodules directly (e.g., package.signals)
            if module.__name__.split('.')[-1] == name:
                obj = module
            else:
                raise

        # Cache it in globals to optimize subsequent lookups
        globals()[name] = obj
        return obj

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# PEP 562: Enables tab-completion in Jupyter/IPython and dir() support
def __dir__():
    return sorted(list(globals().keys()) + __all__)
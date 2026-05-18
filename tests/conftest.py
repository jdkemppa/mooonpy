# -*- coding: utf-8 -*-
"""
Pin the test session to the mooonpy in THIS worktree.

`mooonpy` is installed/importable from the main repo checkout, so without
this a full ``pytest tests/`` run would import that copy as soon as any test
module does ``import mooonpy`` -- and the worktree changes under test would
never be exercised. conftest.py is imported by pytest before any test module,
so prepending the worktree ``src`` here (and evicting any pre-imported
mooonpy) guarantees every test runs against the worktree code.
"""

import os
import sys

_WORKTREE_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))

if sys.path[0] != _WORKTREE_SRC:
    sys.path.insert(0, _WORKTREE_SRC)

# Drop any mooonpy already imported from elsewhere so the worktree copy wins.
for _name in [m for m in list(sys.modules)
              if m == "mooonpy" or m.startswith("mooonpy.")]:
    del sys.modules[_name]

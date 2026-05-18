# -*- coding: utf-8 -*-
"""
Fast, fully independent copy of a Molspace and its containers.

Lazily imported by the ``.copy()`` methods (the default path). The result is
fully independent of the source -- see ``tests/test_molspace_copy.py`` for
the independence/correctness matrix. ``copy.deepcopy`` gives the same result
but is ~6-8x slower (generic per-object walk over ~10^5 slotted objects), so
it is only a cross-check fallback.

Frozen by design, shared by reference (NOT duplicated):
  * the generated value classes (Atom/Bond/Angle/Dihedral/Improper) -- their
    ``__slots__`` are fixed at construction, so every copy reuses the SAME
    class object;
  * read-only lookup config: ``atoms.styles`` (factory + defaults) and
    ``ptable``.

Everything holding per-system mutable state (containers, value objects,
mutable attrs, Box, clusters) is rebuilt fresh, and ``clusters._molspace``
is re-pointed at the new Molspace.
"""

from collections import defaultdict
from copy import deepcopy

from mooonpy.molspace.clusters import Cluster, Clusters
from mooonpy.molspace.force_field import Coefficients, ForceField, Parameters

_IMMUTABLE = (int, float, str, bool, bytes, complex, type(None))


def _dup(v):
    """Independent copy of an attribute value (fast path for scalars)."""
    t = type(v)
    if t in _IMMUTABLE:
        return v
    if isinstance(v, type):
        return v  # share class objects (e.g. generated factory classes)
    if t is list:
        return [_dup(x) for x in v]
    if t is tuple:
        return tuple(_dup(x) for x in v)
    if t is set:
        return {_dup(x) for x in v}
    if t is frozenset:
        return frozenset(_dup(x) for x in v)
    if t is dict:
        return {k: _dup(x) for k, x in v.items()}
    if t is defaultdict:
        new = defaultdict(v.default_factory)
        for k, x in v.items():
            new[k] = _dup(x)
        return new
    return deepcopy(v)  # numpy arrays / unknown objects: correct over fast


def _copy_slotted(obj):
    """Rebuild a __slots__ value object, reusing its (fixed) class."""
    cls = type(obj)
    new = cls.__new__(cls)
    for s in cls.__slots__:
        setattr(new, s, _dup(getattr(obj, s)))
    return new


def _copy_value_dict(src, share_attrs=()):
    """
    Rebuild a dict-subclass container (Atoms/Bonds/Angles/Dihedrals/Impropers)
    of __slots__ value objects. Instance attrs are duplicated, except names in
    `share_attrs` which are shared by reference (read-only config).
    """
    cls = type(src)
    new = cls.__new__(cls)
    for k, v in src.__dict__.items():
        setattr(new, k, v if k in share_attrs else _dup(v))
    for k, obj in src.items():
        dict.__setitem__(new, k, _copy_slotted(obj))
    return new


def _copy_box(src):
    new = type(src).__new__(type(src))
    new.__dict__.update({k: _dup(v) for k, v in src.__dict__.items()})
    return new


def _copy_ff(src):
    new = ForceField.__new__(ForceField)
    for k, v in src.__dict__.items():
        if isinstance(v, Coefficients):
            nc = Coefficients(v.keyword)
            for ak, av in v.__dict__.items():
                setattr(nc, ak, _dup(av))
            for t, param in v.items():
                np_ = Parameters.__new__(Parameters)
                np_.__dict__.update(
                    {ak: _dup(av) for ak, av in param.__dict__.items()}
                )
                dict.__setitem__(nc, t, np_)
            setattr(new, k, nc)
        else:
            setattr(new, k, _dup(v))
    return new


def _copy_clusters(src, new_molspace):
    new = Clusters.__new__(Clusters)
    for k, v in src.__dict__.items():
        # _molspace must point at the NEW Molspace, never the old one
        setattr(new, k, new_molspace if k == '_molspace' else _dup(v))
    for molid, cl in src.items():
        nc = Cluster.__new__(Cluster)
        for s in Cluster.__slots__:
            setattr(nc, s, _dup(getattr(cl, s)))
        dict.__setitem__(new, molid, nc)
    return new


def fast_copy(molspace):
    """Return a fast, fully independent deep copy of ``molspace``."""
    new = molspace.__class__.__new__(molspace.__class__)

    for k, v in molspace.__dict__.items():
        if k == 'atoms':
            na = _copy_value_dict(v, share_attrs=('styles',))
            na.box = _copy_box(v.box)
            new.atoms = na
        elif k in ('bonds', 'angles', 'dihedrals', 'impropers'):
            setattr(new, k, _copy_value_dict(v))
        elif k == 'ff':
            new.ff = _copy_ff(v)
        elif k == 'clusters':
            continue  # built last; needs the new Molspace reference
        elif k == 'ptable':
            new.ptable = v  # static reference data: share
        else:
            setattr(new, k, _dup(v))

    new.clusters = _copy_clusters(molspace.clusters, new)
    return new

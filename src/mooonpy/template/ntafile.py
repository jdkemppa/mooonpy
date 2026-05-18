# -*- coding: utf-8 -*-
"""
NTA-file writer for atom_typing/all2lmp consumption.

Format::

    mooonpy NTA file for <basename> at <date>
    style id

    <atom_id> <nta_string>
    <atom_id> <nta_string>
    ...

The file uses ``style id`` (per-atom-id NTA assignments). Atom NTA strings
come from ``atom.comment`` after canonicalization (only the NTA token is
kept, any verbose ``atom_typing``-style preamble is stripped).
"""
import os
from datetime import datetime


def extract_nta(comment):
    """
    Pull the NTA token out of an atom_typing comment string.

    atom_typing writes verbose comments like
    ``"C        Sp2     avg-angle: 120.000000      type: cp"``. After a
    ``ReactionTemplate.changed_typ`` rewrite the field may already be just
    ``"cp"``. This handles both: returns ``"cp"`` either way, or ``""`` for
    an empty comment.
    """
    if not comment:
        return ''
    parts = comment.split()
    if 'type:' in parts:
        return parts[parts.index('type:') + 1]
    return parts[-1]


def canonicalize_ntas(mol):
    """
    Reduce every ``atom.comment`` in ``mol`` to just its NTA token, in place.

    Idempotent: running it on a Molspace whose comments are already plain
    NTAs is a no-op.
    """
    for atom in mol.atoms.values():
        atom.comment = extract_nta(atom.comment)


def write_ntafile(filename, mol, edge_external_ntas=None):
    """
    Write a per-atom NTA assignment file for ``mol``.

    :param filename: Output path.
    :type filename: str or :class:`Path`
    :param mol: Molspace whose ``atom.comment`` carries the NTA tokens. The
                comment is canonicalized (verbose ``atom_typing`` preambles
                stripped) before writing.
    :param edge_external_ntas: Optional ``{edge_atom_id: [nta, ...]}`` mapping
        each template edge atom to the NTA strings of its bonded neighbors
        *outside* the template. When given (and non-empty), an ``edge id``
        section is appended so downstream charge-equilibration can account
        for the severed bonds. Empty lists are skipped.
    :type edge_external_ntas: dict[int, list[str]] or None
    """
    basename = os.path.basename(str(filename))
    stamp = datetime.now().strftime('%d-%b-%y')
    with open(filename, 'w') as f:
        f.write(f'mooonpy NTA file for {basename} at {stamp}\n')
        f.write('style id\n')
        f.write('\n')
        for id_ in sorted(mol.atoms):
            atom = mol.atoms[id_]
            nta = extract_nta(atom.comment)
            f.write(f'{id_} {nta}\n')

        if edge_external_ntas:
            rows = {eid: ntas for eid, ntas in edge_external_ntas.items() if ntas}
            if rows:
                f.write('\nedge id\n')
                for eid in sorted(rows):
                    f.write('{:5d}\t'.format(eid))
                    for nn in rows[eid]:
                        f.write('{:3} '.format(nn))
                    f.write('\n')

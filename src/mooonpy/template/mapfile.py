# -*- coding: utf-8 -*-
"""
Reaction map-file writer for LAMMPS ``fix bond/react``.
"""
import os
from datetime import datetime


def write_mapfile(filename, rxn_map):
    """
    Write a reaction map file.

    :param filename: Output path.
    :type filename: str or Path
    :param rxn_map: Map data structure. Required keys:

        * ``'Comment'`` (str): free-form comment for the header.
        * ``'InitiatorIDs'`` (list[int]): the two initiator atom ids.
        * ``'EdgeIDs'`` (list[int]): may be empty.
        * ``'DeleteIDs'`` (list[int]): may be empty.
        * ``'Equivalences'`` (list[tuple[int, int]]): ``(pre_id, post_id)`` pairs.
        * ``'CustomChargesQedge'`` (list[int], optional): if non-empty, emits a
          ``CustomChargesQedge`` section so ``fix bond/react`` knows which
          atoms get charges remapped.

    :type rxn_map: dict
    """
    comment = rxn_map.get('Comment', '')
    edges = rxn_map.get('EdgeIDs', [])
    deletes = rxn_map.get('DeleteIDs', [])
    equivs = rxn_map['Equivalences']
    initiators = rxn_map['InitiatorIDs']
    qedge = rxn_map.get('CustomChargesQedge', [])

    basename = os.path.basename(str(filename))

    with open(filename, 'w') as f:
        stamp = datetime.now().strftime('%d-%b-%y')
        f.write(f'mooonpy map file for {basename} at {stamp} {comment}\n')
        f.write('\n')

        if len(edges) > 0:
            f.write(f'{len(edges)} edgeIDs\n')
        f.write(f'{len(equivs)} equivalences\n')
        if len(deletes) > 0:
            f.write(f'{len(deletes)} deleteIDs\n')
        if len(qedge) > 0:
            f.write(f'{len(qedge)} customcharges\n')

        f.write('\n')

        f.write('InitiatorIDs\n')
        f.write('\n')
        for id_ in initiators:
            f.write(f'{id_}\n')
        f.write('\n')

        if len(edges) > 0:
            f.write('EdgeIDs\n')
            f.write('\n')
            for id_ in edges:
                f.write(f'{id_}\n')
            f.write('\n')

        if len(deletes) > 0:
            f.write('DeleteIDs\n')
            f.write('\n')
            for id_ in deletes:
                f.write(f'{id_}\n')
            f.write('\n')

        if len(qedge) > 0:
            f.write('CustomChargesQedge\n')
            f.write('\n')
            for id_ in qedge:
                f.write(f'{id_}\n')
            f.write('\n')

        f.write('Equivalences\n')
        f.write('\n')
        for pre_id, post_id in equivs:
            f.write(f'{pre_id} {post_id}\n')

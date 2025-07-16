# -*- coding: utf-8 -*-



def find(mol, periodicity):
    cutoff = 3
    domains, pairs = mol.compute_pairs(cutoff=cutoff, periodicity=periodicity)
    print(pairs)
    return

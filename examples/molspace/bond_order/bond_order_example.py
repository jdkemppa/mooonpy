# -*- coding: utf-8 -*-
import mooonpy as mp
from mooonpy.tools.file_utils import Path, smart_open
from mooonpy.molspace.force_field import ReaxFF

test_ff = Path('..\\..\\..\\..\\random\\ffield_CHON_for_PAN_PBO.reax')
if not test_ff:
    raise FileNotFoundError
reax = ReaxFF(test_ff)  # read file

mp.rcParams['molspace.astyles'] = ['full', 'reax']

mol = mp.Molspace(Path('..\\WRITE.data'))  # single DETDA
mol.update_elements()

domains, pairs = mol.compute_pairs(3, periodicity='xxx')  # non periodic test case

print('IDs, Types, BO')
for key, pair in pairs.items():
    atom1 = mol.atoms[key[0]]
    atom2 = mol.atoms[key[1]]

    BO_prime = reax.pair_bo_prime(atom1.element, atom2.element, pair.distance)
    bond_key = tuple(reax.sort_ele([atom1.element, atom2.element]))
    # ^ sorted in CHON order based on file ('C', 'H')
    BO_prime = sum(BO_prime)  # 3 parts needed for later correction
    pair.BO_prime = BO_prime
    atom1.BO_total += BO_prime
    atom2.BO_total += BO_prime
    if BO_prime < 0.3: continue
    print(key, bond_key, BO_prime)

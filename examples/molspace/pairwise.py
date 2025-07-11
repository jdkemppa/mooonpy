# -*- coding: utf-8 -*-
"""
@author: Tristan Muzzy
"""
import mooonpy
import matplotlib.pyplot as plt
import numpy as np


# file = mooonpy.Path('../../../big_files/input/relax_EPON_862_pxld_88.2_PCFF_replicate_1_FF_PCFF.data')


# file = mooonpy.Path('../../../big_files/input/B19a-90_R01_04_anneal-017_PBZ.last.data')
# # file = mooonpy.Path('../../../big_files/input/PBZ_old_merged.data')


file = mooonpy.Path('../../../../../Research/PBZ/anneal19/rep*/B19a-90_R*_04_anneal-017_PBZ.last.data')


type_labels = {1:'c', 2:'c2', 3:'c2', 4:'c3', 5:'cp', 6:'cp', 7:'hc', 8:'ho', 9:'nn',10:'oh',11:'op'} # PBZ
bond_hist_total, angle_hist_total, dihedral_hist_total, improper_hist_total = {},{},{},{}

for topofile in file:

    skip_flag = False
    for rep in ['06','07','08','09','10']:
        if rep in topofile:
            skip_flag = True
    if skip_flag: continue
    print(topofile)

    mol = mooonpy.Molspace(topofile)
    mol.add_type_labels(type_labels) # PBZ had no type label
    # mol.atoms.wrap()

    bond_hist, angle_hist, dihedral_hist, improper_hist = mol.compute_BADI_by_type()

    for type_, hist in bond_hist.items():
        if type_ in bond_hist_total:
            bond_hist_total[type_] += hist
        else:
            bond_hist_total[type_] = hist

    for type_, hist in angle_hist.items():
        if type_ in angle_hist_total:
            angle_hist_total[type_] += hist
        else:
            angle_hist_total[type_] = hist

    for type_, hist in dihedral_hist.items():
        if type_ in dihedral_hist_total:
            dihedral_hist_total[type_] += hist
        else:
            dihedral_hist_total[type_] = hist

    for type_, hist in improper_hist.items():
        if type_ in improper_hist_total:
            improper_hist_total[type_] += hist
        else:
            improper_hist_total[type_] = hist


fig, axs = plt.subplots(2,2)

bin_scale = 10

for type_, hist in bond_hist_total.items():
    bins = round(len(hist)/bin_scale)
    if bins <1: continue
    # if 'c2' not in type_: continue
    if 'h'  in type_: continue
    axs[0,0].hist(hist,bins = bins, histtype='step', label=type_)
axs[0,0].legend()
axs[0,0].grid()
axs[0,0].set(title='Bonds')

for type_, hist in angle_hist_total.items():
    bins = round(len(hist)/bin_scale)
    if bins <1: continue
    if 'h' in type_: continue
    axs[0,1].hist(hist,bins = bins, histtype='step', label=type_)
axs[0,1].legend()
axs[0,1].grid()
axs[0,1].set(title='Angles')

for type_, hist in dihedral_hist_total.items():
    bins = round(len(hist)/bin_scale)
    if bins <1: continue
    types = type_.split('-')
    if 'hc' in types: continue
    if 'ho' in types: continue
    if 'c3' in types: continue
    if 'c' in types: continue
    if 'op' in types: continue
    axs[1,0].hist(hist,bins = bins, histtype='step', label=type_)
axs[1,0].legend(ncol=2)
axs[1,0].grid()
axs[1,0].set(title='Dihedrals')


for type_, hist in improper_hist_total.items():
    bins = round(len(hist)/bin_scale)
    if bins <1: continue
    types = type_.split('-')
    axs[1,1].hist(hist,bins = bins, histtype='step', label=type_)
axs[1,1].legend(ncol=2)
axs[1,1].grid()
axs[1,1].set(title='Impropers')

# -*- coding: utf-8 -*-
"""
Plot BADI distributions
@author: Tristan Muzzy
"""
import mooonpy
import matplotlib.pyplot as plt
import numpy as np

# file = mooonpy.Path('../../../big_files/input/relax_EPON_862_pxld_88.2_PCFF_replicate_1_FF_PCFF.data')
# file = mooonpy.Path('../../../big_files/input/B19a-90_R01_04_anneal-017_PBZ.last.data')
# # file = mooonpy.Path('../../../big_files/input/PBZ_old_merged.data')

file = mooonpy.Path('../../../../../Research/PBZ/anneal19/rep*/B19a-90_R*_04_anneal-017_PBZ.last.data')
type_labels = {1: 'c', 2: 'c2', 3: 'c2', 4: 'c3', 5: 'cp', 6: 'cp', 7: 'hc', 8: 'ho', 9: 'nn', 10: 'oh',
               11: 'op'}  # PBZ

# %% Loop through matches
bond_hist_total, angle_hist_total, dihedral_hist_total, improper_hist_total = {}, {}, {}, {}
for topofile in file:
    ## Ignore small reps
    skip_flag = False
    for rep in ['06', '07', '08', '09', '10']:
        if rep in topofile:
            skip_flag = True
    if skip_flag: continue

    print(topofile)
    mol = mooonpy.Molspace(topofile)
    mol.add_type_labels(type_labels)  # PBZ had no type label section
    # mol.atoms.wrap()

    bond_hist, angle_hist, dihedral_hist, improper_hist = mol.compute_BADI_by_type('ppp', comp_improper=False)

    ## combine histos
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

# %% Plots
fig, axs = plt.subplots(2, 2)

bin_scale = 10

for type_, hist in bond_hist_total.items():
    bins = round(len(hist) / bin_scale)
    if bins < 1: continue
    # if 'c2' not in type_: continue
    # if 'h'  in type_: continue
    label = type_
    if mol.ff.has_type_labels:
        for type_, coeff in mol.ff.bond_coeffs.items():
            if label == coeff.type_label:
                label = '{:<7} {:>3.2f}Å'.format(label, coeff.coeffs[0])
                break
    axs[0, 0].hist(hist, bins=bins, histtype='step', label=label)
axs[0, 0].legend()
axs[0, 0].grid()
axs[0, 0].set(title='Bonds')

for type_, hist in angle_hist_total.items():
    bins = round(len(hist) / bin_scale)
    if bins < 1: continue
    if 'h' in type_: continue

    label = type_
    if mol.ff.has_type_labels:
        for type_, coeff in mol.ff.angle_coeffs.items():
            if label == coeff.type_label:
                label = '{:<10} {:>5.2f}°'.format(label, coeff.coeffs[0])
                break
    axs[0, 1].hist(hist, bins=bins, histtype='step', label=label)
axs[0, 1].legend()
axs[0, 1].grid()
axs[0, 1].set(title='Angles')

for type_, hist in dihedral_hist_total.items():
    bins = round(len(hist) / bin_scale)
    if bins < 1: continue
    types = type_.split('-')
    if 'hc' in types: continue
    if 'ho' in types: continue
    if 'c3' in types: continue
    if 'c' in types: continue
    if 'op' in types: continue

    label = type_
    if mol.ff.has_type_labels:
        for type_, coeff in mol.ff.dihedral_coeffs.items():
            if label == coeff.type_label:
                if min(coeff.coeffs) == 0 and max(coeff.coeffs) == 0:
                    label = '{:>13} N/A'.format(label)
                else:
                    minima = []
                    if coeff.coeffs[0] > 0:
                        minima.append(0 + coeff.coeffs[1])
                    elif coeff.coeffs[0] < 0:
                        minima.append(180 + coeff.coeffs[1])

                    if coeff.coeffs[2] > 0:
                        minima += [0 + coeff.coeffs[3], 180 + coeff.coeffs[3]]
                    elif coeff.coeffs[2] < 0:
                        minima += [-90 + coeff.coeffs[3], 90 + coeff.coeffs[3]]

                    if coeff.coeffs[4] > 0:
                        minima += [-120 + coeff.coeffs[5], 0 + coeff.coeffs[5], 120 + coeff.coeffs[5]]
                    elif coeff.coeffs[4] < 0:
                        minima += [-60 + coeff.coeffs[5], 60 + coeff.coeffs[5], 180 + coeff.coeffs[5]]

                    label = '{:<13} {:}°'.format(label, minima)
                break
    axs[1, 0].hist(hist, bins=bins, histtype='step', label=label)
axs[1, 0].legend(ncol=2)
axs[1, 0].grid()
axs[1, 0].set(title='Dihedrals', xlim=[-180, 180])

for type_, hist in improper_hist_total.items():
    bins = round(len(hist) / bin_scale)
    if bins < 1: continue
    types = type_.split('-')

    label = type_
    if mol.ff.has_type_labels:
        for type_, coeff in mol.ff.improper_coeffs.items():
            if label == coeff.type_label:
                label = '{:<13} {:<5.2f}'.format(label, coeff.coeffs[0])
                break
    axs[1, 1].hist(hist, bins=bins, histtype='step', label=label)
axs[1, 1].legend(ncol=2)
axs[1, 1].grid()
axs[1, 1].set(title='Impropers')

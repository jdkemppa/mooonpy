# -*- coding: utf-8 -*-
"""
This is an example of using mooonpy for inserting Li-ions
into a graphite model. You will find details on:
    - How to add FF parameters such as Masses and Pair Coeffs
    - How to find rings in a molecular system
    - How to add atoms to a molecular system
    - How re-wrap all atoms of the newly generated molecular system
"""
import mooonpy as mp
import numpy as np
import random
import os


##########
# Inputs #
##########
# LAMMPS data file to read
file = '../EPON_862/relax_PCFF_IFF_Pitch_rep_1_Graphite_AB.data'

# Lithiation value (0 < x < 1; where 0 is 0% lithiation and 1 is 
# 100% lithiation). A system that is 100% lithiated has a Li-ion
# place in every thing center
lithiation = 1 #0.5

# Seed for random lithiuim ion placement (0 seed will use system time)
seed = 12345

# Distance above ring to place Li-ion (~half the expected d002 spacing)
distance = 3.354/2




####################################
# Start adding to Li-ions to model #
####################################
# Step1: Generate a Molspace instance and add FF parameters
graphite = mp.Molspace(filename=file, astyles=['full', 'charge'])

# Find maxmimum atomTypeID and add Li mass and pair coeffs
li_atom_type = max(graphite.ff.masses.keys()) + 1

# Add to masses 
params = graphite.ff.coeffs_factory()
params.comment = 'li+'
params.type_label = 'li+'
params.coeffs = [6.941]
graphite.ff.masses[li_atom_type] = params


# Add to pair coeffs
params = graphite.ff.coeffs_factory()
params.style = 'lj/class2/coul/long'
params.comment = 'li+'
params.type_label = 'li+'
params.coeffs = [0.000544331, 2.87]
graphite.ff.pair_coeffs[li_atom_type] = params


# Step 2: Find avaiable locations to insert Li-ions (center of 6-member rings)
rings = graphite.find_rings(ring_sizes=(6,))
print('Number of 6-member rings found: {}'.format(len(rings)))

lx, ly, lz = graphite.atoms.box.get_lengths()
print('Simulation cell lengths Lx, Ly, Lz: ', lx, ly, lz)

# set max_x, max_y, max_z w/ minimum image convention
max_x = lx/2; max_y = ly/2; max_z = lz/2
locations = [] # [(x, y, z), ... nlocations]
for n, ring in enumerate(rings):
    if len(ring) != 6: continue
    anchor_atom = min(ring)
    anchor_index = ring.index(anchor_atom)
    atom = graphite.atoms[anchor_atom]
    anchor_x = atom.x
    anchor_y = atom.y
    anchor_z = atom.z
    
    # Unwrap atoms
    xs = [anchor_x]
    ys = [anchor_y]
    zs = [anchor_z]
    for atom_id in ring:
        if atom_id == anchor_atom: continue
        atom = graphite.atoms[atom_id]
        x = atom.x
        y = atom.y
        z = atom.z

        diffx = anchor_x - x
        diffy = anchor_y - y
        diffz = anchor_z - z
        
        # Apply minimum image convention
        if abs(diffx) > max_x:
            x += np.sign(diffx)*lx
        if abs(diffy) > max_y:
            y += np.sign(diffy)*ly
        if abs(diffz) > max_z:
            z += np.sign(diffz)*lz
        
        # Save atom positions
        xs.append(x)
        ys.append(y)
        zs.append(z)
        
    # Find center of atoms and shift in Z-direction,
    # where future Li-ion's should be place
    xc = sum(xs)/len(xs)
    yc = sum(ys)/len(ys)
    zc = sum(zs)/len(zs) + distance
    locations.append((xc, yc, zc))
    
    
# Step 3: Insert Li-ion's up until the maximum lithiation value
if seed > 0: random.seed(seed)
if lithiation < 1.0:
    n_ions = int(np.floor(lithiation*len(locations)))
else:
    n_ions = len(locations)

atom_id = max(graphite.atoms.keys())
n_inserted = 0
while n_inserted < n_ions:
    random_index = random.randint(0, len(locations)-1)
    x, y, z = locations[random_index]
    del locations[random_index]
    
    atom_id += 1
    atom = graphite.atoms.styles.atom_factory()
    atom.x = x
    atom.y = y
    atom.z = z
    atom.q = 1.0
    atom.id = atom_id
    atom.type = li_atom_type
    atom.comment = 'li+'
    graphite.atoms[atom_id] = atom
    n_inserted += 1

print('Inserted {} Li-ions of {} requested'.format(n_inserted, n_ions))

        
# Step 4: Wrap atoms and write LAMMPS datafile
basename = os.path.splitext(file)[0]
graphite.atoms.wrap()
graphite.write_files('{}_Li-ions.data'.format(basename), atom_style='full')
# -*- coding: utf-8 -*-



def find(mol, periodicity):
    # Find all elements in system
    elements = set([])
    for atom_type, mass in mol.ff.masses.items():
        element = mass.element
        if element == '':
            print(f'\nERROR {__file__}')
            raise Exception(f'Mass {atom_type} does not have element set. Update elements via "my_molecule.update_elements()"')
        elements.add(element)
    
    # Find per-element information
    valences = {} # {'element': valence}
    atom_sizes = {} # {'element': atom_size}
    for element in elements:
        ptable = mol.ptable.elements[element]
        valences[element] = ptable.valence
        atom_sizes[element] = ptable.radii['vdw']
    
    # Set cutoff as twice that largest atoms vdw radii
    cutoff = 2*atom_sizes[max(atom_sizes, key=atom_sizes.get)]
    
    # Find pairs of atoms within cutoff and generate graph
    domains, pairs = mol.compute_pairs(cutoff=cutoff, periodicity=periodicity)
    graph = {i:[] for i in mol.atoms}
    for (id1, id2), info in pairs.items():    
        graph[id1].append(id2)
        graph[id2].append(id1)
    
    # Start cutting down pairs
    bonds, flagged_bonds = set([]), set([])
    for id1 in graph:       
        atom1 = mol.atoms[id1]
        radii1 = atom_sizes[atom1.element]
        max_nb = valences[atom1.element]

        # Creating dictionary to sort in decending order
        distances = {}
        for id2 in graph[id1]:
            atom2 = mol.atoms[id2]
            radii2 = atom_sizes[atom2.element]
            cutoff = radii1 + radii2
            bond = tuple(sorted([id1, id2]))
            distance = pairs[bond].distance
            if distance <= cutoff: 
                distances[bond] = distance
        distances = dict( sorted(distances.items(), key=lambda x:x[1]) )
        
        # Creating bonds based on inputs values
        count_nb = 0;
        for bond in distances:
            count_nb += 1
            # If number of bonds is less than specified
            if count_nb <= max_nb:
                bonds.add(bond)
            # If bond does not meet criteria flag the bond to not create the bond
            else:
                flagged_bonds.add(bond)
                
    # Removing any flagged bonds if they were created      
    bonds = list(set(bonds) - set(flagged_bonds)) 
   
    # Sorting bonds and flagged bonds                
    bonds = sorted(bonds)
    flagged_bonds = sorted(flagged_bonds) 
    
    # Clear any existing bonds and add new bonds to system
    mol.bonds.clear()
    for (id1, id2) in bonds:
        element1 = mol.atoms[id1].element
        element2 = mol.atoms[id2].element
        distance = pairs[(id1, id2)].distance
        comment = '{}-{}: {}'.format(element1, element2, distance)
        
        bond = mol.bonds.bond_factory()
        bond.ordered = [id1, id2]
        bond.comment = comment
        bond.type = 1
        mol.bonds.id_set(bond, id1, id2)

    return

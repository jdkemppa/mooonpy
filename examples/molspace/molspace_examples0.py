# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 16:03:17 2025

@author: jdkem
"""
import mooonpy
#from mooonpy import guis, tools, molspace

#import mooonpy
#from mooonpy import molspace as ms
#from mooonpy import thermospace as ts

# print(mooonpy.molspace.molspace.Molspace())
# print(mooonpy.Molspace())

# print(mooonpy.molspace.doc_examples.add(1,3))
# print(mooonpy.DocExamples.add(1,3))


#print('full : ', mooonpy) #.doc_examples.add(1, 2))
#print('alias: ', ms.hw.hello_world())


#print('full : ', mooonpy.thermospace.multiply(x=5, y=10))
#print('alias: ',  ts.lw.multiply(x=5, y=10))


# file = '../EPON_862/all2lmp_Outputs/detda_typed_IFF.data'
# #file = '../EPON_862/Graphite_AB_relaxed.data'
file = '../EPON_862/dgebf_typed_IFF_merged.data'
file = '../EPON_862/pyrene_typed_IFF.data'
file = '../EPON_862/ortho_catechol_mxene_sheet_IFF.data'

# #file = '../EPON_862/Cellulose-supercell_morse_IFF.data'
# file = '../EPON_862/system1_cell_replicate.data'
#file = '../EPON_862/detda_typed_IFF_merged.data'
file = '../EPON_862/detda_typed_IFF_hybrid_class2_morse.data'
#file = '../EPON_862/Graphite_AB_relaxed.data'


mooonpy.rcParams['color'] = 'green'


molecule = mooonpy.Molspace(filename=file)#, astyles=['full', 'charge'], dsect=['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities'])
molecule.write_files('WRITE.data', atom_style='full')
molecule.write_files('WRITE.ff.script')

    
rings = molecule.find_rings(ring_sizes=(3,4,5,6,7))
print(rings)

molecule.update_per_atom_element()
molecule.bonds_from_distances(periodicity='fff')

print('\n\nMasses')
for i in molecule.ff.masses:
    mass = molecule.ff.masses[i]
    print(i, mass.comment, mass.element)

print('\n\nAtoms')
for i in molecule.atoms:
    atom = molecule.atoms[i]
    print(i, atom.type, atom.comment, atom.element)
    

    

    
    

    

        
        

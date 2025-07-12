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
#file = '../EPON_862/ortho_catechol_mxene_sheet_IFF.data'

# #file = '../EPON_862/Cellulose-supercell_morse_IFF.data'
# file = '../EPON_862/system1_cell_replicate.data'
#file = '../EPON_862/detda_typed_IFF_merged.data'
file = '../EPON_862/detda_typed_IFF_hybrid_class2_morse.data'
#file = '../EPON_862/Graphite_AB_relaxed.data'


mooonpy.rcParams['color'] = 'green'


molecule = mooonpy.Molspace(filename=file)#, astyles=['full', 'charge'], dsect=['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities'])
molecule.write_files('WRITE.data', atom_style='full')
molecule.write_files('WRITE.ff.script')

#lines = molecule.ff.get_per_line_styles('bond_coeffs')
#print(lines)



graph = mooonpy.molspace.graph_theory.interface.generate_graph(molecule)
#print(graph)

# rings = mooonpy.molspace.graph_theory.graph.find_rings(graph, ring_sizes=(3,4,5,6,7))
# print(len(rings))

def call():
    rings = mooonpy.molspace.graph_theory.ring_analysis.find_rings(graph, ring_sizes=(3,4,5,6,7))
    print()
    print(rings)
    print(len(rings))
    
    rings = molecule.find_rings(ring_sizes=(3,4,5,6,7))
    print(rings)
    
import timeit
time = timeit.timeit(stmt=call, number=1)
print(time)



if __name__ == '__main__':
    import timeit
    
    # file = 'EPON_862/detda_typed_IFF_merged.data'
    # file = 'EPON_862/system1_cell_replicate.data'
    # def call_mooonpy():
    #     m = mooonpy.Molspace(filename=file, read='mooonpy', astyles=['all', 'full'])
    #     m.write_files('WRITE.data', atom_style='full')
    
    # number = 1
    # print('\n\n')
    # mooonpy_time = timeit.timeit(stmt=call_mooonpy, number=number)
    # print(f'mooonpy read time  : {mooonpy_time} seconds for {number} runs on 100,000 atom system')
    
    

    
    

    

        
        

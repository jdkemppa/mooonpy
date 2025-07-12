# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 10:40:21 2025

@author: jdkem

https://ptable.com/?lang=en#Properties
"""

class Element:
    def __init__(self):
        # mass are in AMU's
        self.masses: list[float] = [0.0]
        
        # Values are in angstrom's
        self.radii: dict[float] = {'calculated': 0.0,
                                   'empirical':  0.0,
                                   'covalent':   0.0,
                                   'vdw':        0.0,
                                   'ff.ReaxFF':  0.0,
                                   'ff.REBO':    0.0}
        
        # Might not really care about this info
        self.elevels: str = ''
        

class Elements:
    def __init__(self):
        self.elements = {} # {'element': Element-object}
        
        carbon = Element()
        carbon.masses = [12.011, 10.01115] # 10.01115 is for IFF's cg1/cge atom types
        carbon.radii = {'calculated': 0.67,
                        'empirical':  0.70,
                        'covalent':   0.77,
                        'vdw':        1.70}
        self.elements['C'] = carbon
        
        
        hydrogen = Element()
        hydrogen.masses = [1.008, 1.0] # 1.0 is for IFF's cg1/cge atom types
        hydrogen.radii = {'calculated': 0.53,
                          'empirical':  0.25,
                          'covalent':   0.37,
                          'vdw':        1.20}
        self.elements['H'] = hydrogen
        
    def mass2element(self, mass):
        mass_diffs = {} # {'element':minimum-difference in masses}
        for elem in self.elements:
            masses = self.elements[elem].masses
            diffs = [abs(mass - elem_mass) for elem_mass in masses]
            mass_diffs[elem] = min(diffs)
        return min(mass_diffs, key=mass_diffs.get)
    
    def element2mass(self, element):
        return self.elements[element].masses[0]
    
    def element2radii(self, element, method='vdw'):
        return self.elements[element].radii[method]

if __name__ == "__main__":
    pt = Elements()
    carbon = pt.elements['C']
    print(carbon.masses, carbon.radii)
    
    print('\n\nMapping mass to element')
    print(pt.mass2element(12))
    print(pt.element2mass('H'))
    print(pt.element2radii('C'))
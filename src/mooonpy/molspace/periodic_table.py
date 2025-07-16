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
        
        # Valence number 
        self.valence: int = 0
        
        # Might not really care about this info
        self.elevels: str = ''
        
# https://ptable.com/?lang=en#Properties
class Elements:
    def __init__(self):
        self.elements = {} # {'element': Element-object}
        
        carbon = Element()
        carbon.masses = [12.011] # Firt index will be the formal mass used when getting a mass per element
        carbon.masses.extend([12.01115, 10.01115])  # Other FF masses (10.01115 is for IFF's cg1/cge atom types)
        carbon.radii = {'calculated': 0.67,
                        'empirical':  0.70,
                        'covalent':   0.77,
                        'vdw':        1.70}
        carbon.valence = 4
        self.elements['C'] = carbon
        
        
        hydrogen = Element()
        hydrogen.masses = [1.008] # Firt index will be the formal mass used when getting a mass per element
        hydrogen.masses.extend([1.0, 1.00782, 1.00797, 1.008, 2.014])  # Other FF masses (1.0 is for IFF's cg1/cge atom types)
        hydrogen.radii = {'calculated': 0.53,
                          'empirical':  0.25,
                          'covalent':   0.37,
                          'vdw':        1.20}
        hydrogen.valence = 1
        self.elements['H'] = hydrogen
        
        oxygen = Element()
        oxygen.masses = [15.999] # Firt index will be the formal mass used when getting a mass per element
        oxygen.masses.extend([14.9994, 15.99491, 15.9994, 16.0])  # Other FF masses
        oxygen.radii = {'calculated': 0.48,
                        'empirical':  0.60,
                        'covalent':   0.73,
                        'vdw':        1.52}
        oxygen.valence = 2
        self.elements['O'] = oxygen
        
        nitrogen = Element()
        nitrogen.masses = [14.007] # Firt index will be the formal mass used when getting a mass per element
        nitrogen.masses.extend([14.0067, 14.00674, 14.01, 14.0])  # Other FF masses
        nitrogen.radii = {'calculated': 0.56,
                          'empirical':  0.65,
                          'covalent':   0.75,
                          'vdw':        1.55}
        nitrogen.valence = 4
        self.elements['N'] = nitrogen
        
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
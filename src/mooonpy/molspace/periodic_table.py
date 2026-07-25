# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 10:40:21 2025

@author: jdkem, tkwav

https://ptable.com/?lang=en#Properties
"""

from typing import Union
from dataclasses import dataclass, field, fields

@dataclass
class Element:
        name: str = field(default_factory=str)
        symbol: str = field(default_factory=str)
        number: int = field(default_factory=int) # atomic number
        masses: list[float] = field(default_factory=list)
        radius: dict[str, float] = field(default_factory=dict) # list of radii types in Angstrom
        valence: float = field(default_factory=float)

        @property
        def mass(self) -> Union[float, None]:
            if self.masses:
                return self.masses[0]
            else:
                return None

        def __repr__(self) -> str:
            return self.symbol


@dataclass
class PeriodicTable:
    elements: dict[str, Element] = field(default_factory=dict)
    _by_number: dict[int, Element] = field(default_factory=dict, init=False, repr=False)

    # Fill in most commonly used elements. Can add extras using PeriodicTable.add_element()
    def __post_init__(self):
        _defaults = [
            Element(name='Carbon',
                    symbol="C",
                    number=6,
                    masses=[12.011, 12.01115, 10.01115],
                    radius={'calculated': 0.67,
                            'empirical': 0.70,
                            'covalent': 0.77,
                            'vdw': 1.70},
                    valence=4),

            Element(name='Hydrogen',
                    symbol="H",
                    number=1,
                    masses=[1.008, 1.0, 1.00782, 1.00797, 1.008, 2.014],
                    radius={'calculated': 0.53,
                            'empirical': 0.25,
                            'covalent': 0.37,
                            'vdw': 1.20},
                    valence=1),

            Element(name='Oxygen',
                    symbol="O",
                    number=8,
                    masses=[15.999, 14.9994, 15.99491, 15.9994, 16.0],
                    radius={'calculated': 0.48,
                            'empirical': 0.60,
                            'covalent': 0.73,
                            'vdw': 1.52},
                    valence=2),

            Element(name='Nitrogen',
                    symbol="N",
                    number=7,
                    masses=[14.007, 14.0067, 14.00674, 14.01, 14.0],
                    radius={'calculated': 0.56,
                            'empirical': 0.65,
                            'covalent': 0.75,
                            'vdw': 1.55},
                    valence=4)
        ]
        for def_element in _defaults:
            self.elements[def_element.symbol] = def_element
            self._by_number[def_element.number] = def_element

    def __str__(self):
        rows = [f"  {e.symbol:<6}{e.name:<12}{e.masses[0]} amu"
                for e in self.elements.values()]
        return f"PeriodicTable Object with ({len(self.elements)} elements):\n" + "\n".join(rows)

    def add_element(self, element: Element):
        self.elements[element.symbol] = element
        self._by_number[element.number] = element

    def get(self, symbol: str) -> Element | None:
        return self.elements.get(symbol)

    from_symbol = get  # aliased

    def from_number(self, number: int) -> Element | None:
        return self._by_number.get(number)

    def to_dict(self):
        return {f.name: getattr(self, f.name) for f in fields(self)}

    def mass2element(self, mass: float):
        # retrieve element of closest mass
        mass_diffs = {}  # {'element':minimum-difference in masses}
        for elem in self.elements:
            masses = self.elements[elem].masses
            diffs = [abs(mass - elem_mass) for elem_mass in masses]
            mass_diffs[elem] = min(diffs)
        return min(mass_diffs, key=mass_diffs.get)

    def element2mass(self, element):
        return self.elements[element].masses[0]

    def element2radius(self, element, method='vdw'):
        return self.elements[element].radius[method]


PERIODIC_TABLE = PeriodicTable()

if __name__ == "__main__":
    pt = PeriodicTable()
    carbon = pt.elements['C']
    print(carbon.masses, carbon.radius)

    print('\n\nMapping mass to element')
    print(pt.mass2element(12))
    print(pt.element2mass('H'))
    print(pt.element2radius('C'))

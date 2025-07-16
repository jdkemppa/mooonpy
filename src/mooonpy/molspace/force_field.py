# -*- coding: utf-8 -*-


class Parameters(object):
    def __init__(self, coeffs):                
        self.coeffs: list = coeffs
        self.style: str = ''
        self.comment: str = ''
        self.type_label: str = ''
        self.element: str = ''


class Coefficients(dict):
    def __init__(self, keyword, *args,  **kwargs):
        super().__init__(*args, **kwargs)
        self.style: str = ''
        self.keyword = keyword
        

class ForceField(object):
    def __init__(self, **kwargs):
        # Build this object with some composition
        self.masses: Coefficients = Coefficients('Masses')           # {type : Parameters-object}
        self.pair_coeffs: Coefficients = Coefficients('Pair Coeffs') # {type : Parameters-object}
        
        self.bond_coeffs: Coefficients = Coefficients('Bond Coeffs')         # {type : Parameters-object}
        self.angle_coeffs: Coefficients = Coefficients('Angle Coeffs')       # {type : Parameters-object}
        self.dihedral_coeffs: Coefficients = Coefficients('Dihedral Coeffs') # {type : Parameters-object}
        self.improper_coeffs: Coefficients = Coefficients('Improper Coeffs') # {type : Parameters-object}

        self.bondbond_coeffs: Coefficients = Coefficients('BondBond Coeffs')                   # {type : Parameters-object}
        self.bondangle_coeffs: Coefficients = Coefficients('BondAngle Coeffs')                 # {type : Parameters-object}
        self.angleangletorsion_coeffs: Coefficients = Coefficients('AngleAngleTorsion Coeffs') # {type : Parameters-object}
        self.endbondtorsion_coeffs: Coefficients = Coefficients('EndBondTorsion Coeffs')       # {type : Parameters-object}
        self.middlebondtorsion_coeffs: Coefficients = Coefficients('MiddleBondTorsion Coeffs') # {type : Parameters-object}
        self.bondbond13_coeffs: Coefficients = Coefficients('BondBond13 Coeffs')               # {type : Parameters-object}        
        self.angletorsion_coeffs: Coefficients = Coefficients('AngleTorsion Coeffs')           # {type : Parameters-object}
        self.angleangle_coeffs: Coefficients = Coefficients('AngleAngle Coeffs')               # {type : Parameters-object}
        
        # Boolean to check if type labels have been read-in
        self.has_type_labels = False
        
        
    def coeffs_factory(self, coeffs=None):
        if coeffs is None: coeffs = []
        return Parameters(coeffs)
    
    def get_per_line_styles(self, coeff):
        lines = {} # {'TypeID':'per-line-style'}
        potential = getattr(self, coeff)
        for i in potential:
            parameters = potential[i]
            lines[i] = parameters.style
        return lines

        


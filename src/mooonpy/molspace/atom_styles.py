# -*- coding: utf-8 -*-
import re
from collections import defaultdict
from copy import deepcopy


def _is_float(string):
    float_re = re.compile(r'^-?\d+(\.\d+)?([eE][-+]?\d+)?$')
    return bool(float_re.match(string))


def _int_str(string):
    if string.isnumeric():
        return int(string)
    else:
        return string

    
def _int_str_float(string):
    if string.isnumeric():
        return int(string)
    elif _is_float(string):
        return float(string)
    else:
        return string


def _make_class(class_name, slots, defaults=None):
    defaults = defaults or {}
    
    class Dynamic:
        __slots__ = slots
        def __init__(self, values=None):
            values = values or {}
            merged = {**defaults, **values}
            for name in self.__slots__:
                setattr(self, name, merged.get(name))
        def copy(self):
            return deepcopy(self)

    Dynamic.__name__ = class_name
    return Dynamic


class Styles:
    # Setup the default values per each attribute
    defaults = defaultdict(str)
    defaults.update({'id': 0,
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0,
                    'q': 0.0,
                    'x0': 0.0,
                    'y0': 0.0,
                    'z0': 0.0,
                    'vx': 0.0,
                    'vy': 0.0,
                    'vz': 0.0,
                    'ix': 0,
                    'iy': 0,
                    'iz': 0,
                    'fx': 0.0,
                    'fy': 0.0,
                    'fz': 0.0,
                    'xu': 0.0,
                    'yu': 0.0,
                    'zu': 0.0,
                    'cv': 0.0,
                    'ed': 0.0,
                    'em': 0.0,
                    'sp': 0.0,
                    'spx': 0.0,
                    'spy': 0.0,
                    'spz': 0.0,
                    'mux': 0.0,
                    'muy': 0.0,
                    'muz': 0.0,
                    'rho': 0.0,
                    'ccN': 0.0,
                    'type': 0,
                    'etag': 0,
                    'esph': 0.0,
                    'area': 0.0,
                    'mass': 0.0,
                    'cs_re': 0.0,
                    'cs_im': 0.0,
                    'molid': 0,
                    'espin': 0,
                    'theta': 0.0,
                    'energy': 0.0,
                    'status': 0,
                    'volume': 0.0,
                    'edpd_cv': 0.0,
                    'density': 0.0,
                    'epsilon': 0.0,
                    'eradius': 0.0,
                    'kradius': 0.0,
                    'cradius': 0.0,
                    'molecule': 0,
                    'bodyflag': 0,
                    'lineflag': 0,
                    'diameter': 0.0,
                    'curvature': 0.0,
                    'edpd_temp': 0.0,
                    'triangleflag': 0,
                    'ellipsoidflag': 0,
                    'template_atom': 0,
                    'template_index': 0,
                    'sub_styleN': 0.0,
                    'comment': '',
                    'element': '',
                    'name': '',
                    'BO_total': 0.0,
                    })

    # Set the function aliases on how to convert a string from a file
    # to a python data type, during the reading of a file.
    func_aliases = defaultdict(lambda: str)
    func_aliases.update({'id': int,
                        'x': float,
                        'y': float,
                        'z': float,
                        'q': float,
                        'x0': float,
                        'y0': float,
                        'z0': float,
                        'vx': float,
                        'vy': float,
                        'vz': float,
                        'ix': int,
                        'iy': int,
                        'iz': int,
                        'fx': float,
                        'fy': float,
                        'fz': float,
                        'xu': float,
                        'yu': float,
                        'zu': float,
                        'cv': float,
                        'ed': float,
                        'em': float,
                        'sp': float,
                        'spx': float,
                        'spy': float,
                        'spz': float,
                        'mux': float,
                        'muy': float,
                        'muz': float,
                        'rho': float,
                        'ccN': float,
                        'type': _int_str,
                        'etag': int,
                        'esph': float,
                        'area': float,
                        'mass': float,
                        'cs_re': float,
                        'cs_im': float,
                        'molid': int,
                        'espin': int,
                        'theta': float,
                        'energy': float,
                        'status': int,
                        'volume': float,
                        'edpd_cv': float,
                        'density': float,
                        'epsilon': float,
                        'eradius': float,
                        'kradius': float,
                        'cradius': float,
                        'molecule': int,
                        'bodyflag': int,
                        'lineflag': int,
                        'diameter': float,
                        'curvature': float,
                        'edpd_temp': float,
                        'triangleflag': int,
                        'ellipsoidflag': int,
                        'template_atom': int,
                        'template_index': int,
                        'sub_styleN': _int_str_float,
                        'comment': str,
                        'element': str,
                        'name': str,
                        'BO_total': float,
                        })

    styles = {}  # {'style' : (attr1, attr2, ...) }      -> (attr1, attr2, ...) sets attr name (e.g. 'id' 'x')
    styles['angle'] = ('id', 'molid', 'type', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['atomic'] = ('id', 'type', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['body'] = ('id', 'type', 'bodyflag', 'mass', 'x', 'y', 'z', 'ix', 'iy', 'iz')

    styles['bond'] = ('id', 'molid', 'type', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['bpm/sphere'] = ('id', 'molid', 'type', 'diameter', 'density', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['charge'] = ('id', 'type', 'q', 'x', 'y', 'z', 'ix', 'iy', 'iz')

    styles['dielectric'] = ('id', 'type', 'q', 'x', 'y', 'z', 'mux', 'muy', 'muz', 'area', 'ed', 'em', 'epsilon',
                            'curvature', 'ix', 'iy', 'iz')
    styles['dipole'] = ('id', 'type', 'q', 'x', 'y', 'z', 'mux', 'muy', 'muz', 'ix', 'iy', 'iz')
    styles['dpd'] = ('id', 'type', 'theta', 'x', 'y', 'z', 'ix', 'iy', 'iz')

    styles['edpd'] = ('id', 'type', 'edpd_temp', 'edpd_cv', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['electron'] = ('id', 'type', 'espin', 'eradius', 'q', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['ellipsoid'] = ('id', 'type', 'ellipsoidflag', 'density', 'x', 'y', 'z', 'ix', 'iy', 'iz')

    styles['full'] = ('id', 'molid', 'type', 'q', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['line'] = ('id', 'molid', 'type', 'lineflag', 'density', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['mdpd'] = ('id', 'type', 'rho', 'x', 'y', 'z', 'ix', 'iy', 'iz')

    styles['molecular'] = ('id', 'molid', 'type', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['peri'] = ('id', 'type', 'volume', 'density', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['rheo'] = ('id', 'type', 'status', 'rho', 'x', 'y', 'z', 'ix', 'iy', 'iz')

    styles['rheo/thermal'] = ('id', 'type', 'status', 'rho', 'energy', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['smd'] = ('id', 'type', 'molecule', 'volume', 'mass', 'kradius', 'cradius', 'x0', 'y0', 'z0', 'x', 'y',
                     'z', 'ix', 'iy', 'iz')
    styles['sph'] = ('id', 'type', 'rho', 'esph', 'cv', 'x', 'y', 'z', 'ix', 'iy', 'iz')

    styles['sphere'] = ('id', 'type', 'diameter', 'density', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['spin'] = ('id', 'type', 'x', 'y', 'z', 'spx', 'spy', 'spz', 'sp', 'ix', 'iy', 'iz')
    styles['tdpd'] = ('id', 'type', 'x', 'y', 'z',
                      'ccN')  # WARNING 'ccN', needs to be expanded to 'cc1', 'cc2', ... so 'tdpd' NOT TRULY SUPPORTED

    styles['template'] = ('id', 'type', 'molid', 'template_index', 'template_atom', 'x', 'y', 'z', 'ix', 'iy',
                          'iz')
    styles['tri'] = ('id', 'molid', 'type', 'triangleflag', 'density', 'x', 'y', 'z', 'ix', 'iy', 'iz')
    styles['wavepacket'] = ('id', 'type', 'q', 'espin', 'eradius', 'etag', 'cs_re', 'cs_im', 'x', 'y', 'z', 'ix',
                            'iy', 'iz')

    styles['hybrid'] = ('id', 'type', 'x', 'y', 'z',
                        'sub_styleN')  # WARNING 'sub_styleN', needs to be expanded to 'sub_style1', 'sub_style2', ... so 'hybrid' NOT TRULY SUPPORTED
    styles['reax'] = ('BO_total',)
    styles['box'] = ('xu', 'yu', 'zu',)
    styles['force'] = ('fx', 'fy', 'fz',)
    styles['custom'] = ()

    styles['_random'] = ('comment', 'element', 'name', 'vx', 'vy', 'vz',)

    # self.styles['_chemistry']     = ('hybrid', 'element', 'rings') # attribute groupings
    ## These are class attributes so they may be used as a lookup externally
    def __init__(self, astyles):
        """
        Set up supported styles for LAMMPS and other file formats (NOTE: for the LAMMPS 'type' attributes they can be an
        int or str due to type labels in a LAMMPS datafile or LAMMPS molecule file - so need to use _int_str function).

        The self.read, self.styles, and self.defaults dictionaries are meant to provide very easy general purpose way
        to extend the ability to read and write LAMMPS styles. However, due to the general purpose nature they are not
        as quick as they could be. Therefore, "hard coded" ways to read and write for speed purposes can be added. These
        will be set as functions with the names:
           def read_STYLE(self, atom, style, data_lst):
           def line_STYLE(self, atom, style):
        where STYLE in def *_STYLE() is a LAMMPS style (e.g. def read_full(), will be the hard coded reader to the
        LAMMPS 'full' atom style). *NOTE: to support any LAMMPS style that style key needs to be in self.read, self.styles,
        and self.defaults dictionaries, even if "hard coded" readers and writers are generated - this is because this
        will tell the code how to generate an instance of the Atom class, with necessary __slots__ and defaults.*
        """

        # Build the dictionaries mapping each key/value pair from above into each style 
        self.read = {}     # {'style' : (int or float or str)}     -> (int or float or str) sets how to read string
        self.defaults = {} # {'style' : (default1, default2, ...)} -> (default1, default2, ...) sets default value for each attr
        for style in self.styles:
            _defaults, _func_aliases = [], []
            for attr in self.styles[style]:
                if attr in Styles.defaults:
                    default = Styles.defaults[attr]
                else:
                    default = 0
                    print(f'\nWARNING {__file__}')
                    print(f'is setting a general default: "{default}" for per-atom attribute: "{attr}"') 
                _defaults.append(default)
                    
                if attr in self.func_aliases:
                    func_alias = self.func_aliases[attr]
                else:
                    func_alias = _int_str_float
                    print(f'\nWARNING {__file__}')
                    print(f'is setting a general function alias for data type conversions: "{func_alias}" for per-atom attribute: "{attr}"') 
                _func_aliases.append(func_alias)
                
            self.defaults[style] = tuple(_defaults)
            self.read[style] = tuple(_func_aliases)
        
        
        # Setup which styles should be built when an atom is generated
        if isinstance(astyles, (tuple, list)):
            build = tuple(astyles)
        else:
            build = tuple([astyles])

        
        # Consolidate all attrs and defaults into a dict and tuple to 
        # be able to initialize an Atom() object as quickly as possible
        self.all_defaults = {} # {'attr-name':default-value}
        #print('Generating atom styles: ')
        for style in self.styles:
            if style in build or 'all' in build or style.startswith('_'):
                #print(' - Building atom style: ', style)
                attrs = self.styles[style]
                defaults = self.defaults[style]
                for attr, default in zip(attrs, defaults):
                    self.all_defaults[attr] = default
        self.all_per_atom = tuple(self.all_defaults.keys())

        
        # Generate an Atom class with necessary slots and defaults, which will
        # be used by atom_factory() method to generate an instance of this class
        class_name = 'Atom'
        slots = self.all_per_atom
        defaults = self.all_defaults
        # print('slots', slots)
        self.Atom = _make_class(class_name, slots, defaults=defaults)

        
    def update_astyles(self, astyles=None):
        # Update defaults
        if astyles is None:
            astyles = self.styles
        for style in astyles:
            attrs = self.styles[style]
            defaults = self.defaults.get(style,Styles.defaults[style])
            print(style, attrs, defaults)
            for attr, default in zip(attrs, defaults):
                self.all_defaults[attr] = default
        self.all_per_atom = tuple(self.all_defaults.keys())
        # print(self.all_per_atom)
        # Re-generate the Atom class with the new defaults
        class_name = 'Atom'
        slots = self.all_per_atom
        defaults = self.all_defaults
        self.Atom = _make_class(class_name, slots, defaults=defaults)
        return
        
    def atom_factory(self):
        return self.Atom()
    
    def atom_fill(self, atom, style, data_lst):
        reader = self.read[style]
        attrs = self.styles[style]
        for string, attr, func in zip(data_lst, attrs, reader):
            value = func(string)
            setattr(atom, attr, value)
        return atom
    
    def read_full(self, atom, style, data_lst):
        atom.id = int(data_lst[0])
        atom.molid = int(data_lst[1])
        atom.type = int(data_lst[2])
        atom.q = float(data_lst[3])
        atom.x = float(data_lst[4])
        atom.y = float(data_lst[5])
        atom.z = float(data_lst[6])
        try:
            atom.ix = int(data_lst[7])
            atom.iy = int(data_lst[8])
            atom.iz = int(data_lst[9])
        except:
            atom.ix = 0
            atom.iy = 0
            atom.iz = 0
        return atom
    
    def read_charge(self, atom, style, data_lst):
        atom.id = int(data_lst[0])
        atom.type = int(data_lst[1])
        atom.q = float(data_lst[2])
        atom.x = float(data_lst[3])
        atom.y = float(data_lst[4])
        atom.z = float(data_lst[5])
        try:
            atom.ix = int(data_lst[6])
            atom.iy = int(data_lst[7])
            atom.iz = int(data_lst[8])
        except:
            atom.ix = 0
            atom.iy = 0
            atom.iz = 0
        return atom
    
    def atom_line(self, atom, style):
        attrs = self.styles[style]
        buffer = 2
        line = ''
        for attr in attrs:
            value = getattr(atom, attr)
            if isinstance(value, (int, str)):
                nchars = len(str(value))
                space = nchars + 2*buffer
                string = '{text:^{s}}'.format(text=value, s=space)
            elif isinstance(value, float):
                value = '{:.16f}'.format(value)
                nchars = len(str(value))
                space = nchars + 2*buffer
                string = '{text:^{s}}'.format(text=value, s=space)
            else:
                string = '{}'.format(value)  
            line += string
        return line
    
    def line_full(self, atom, style):
        pos = '{:^24.16f} {:^24.16f} {:^24.16f}'.format(atom.x, atom.y, atom.z)
        image = '{:^2} {:^2} {:^2}'.format(atom.ix, atom.iy, atom.iz)
        
        line = '{:^6} {:^4} {:^2} {:^20.16f} {} {}'.format(atom.id, atom.molid, atom.type, atom.q, pos, image)
        return line
    
    def line_charge(self, atom, style):
        pos = '{:^24.16f} {:^24.16f} {:^24.16f}'.format(atom.x, atom.y, atom.z)
        image = '{:^2} {:^2} {:^2}'.format(atom.ix, atom.iy, atom.iz)
        
        line = '{:^6} {:^4} {:^20.16f} {} {}'.format(atom.id, atom.type, atom.q, pos, image)
        return line
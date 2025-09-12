# -*- coding: utf-8 -*-
from mooonpy.molspace import _files_io as _files_io
from mooonpy.molspace.atoms import Atoms
from mooonpy.molspace.topology import Bonds, Angles, Dihedrals, Impropers
from mooonpy.molspace.clusters import Clusters
from mooonpy.molspace.force_field import ForceField
from mooonpy.molspace.distance import domain_decomp_13, pairs_from_bonds, pairs_from_domains, ADI_from_bonds, \
    BADI_by_type
from mooonpy.molspace.graph_theory import interface, ring_analysis
from mooonpy.molspace.periodic_table import Elements as Ptable
from mooonpy.molspace.bonds_from_distances import find as find_bonds_from_distances

from mooonpy.tools.file_utils import Path

from mooonpy.rcsetup import rcParams
from copy import deepcopy
import warnings
import os


class Molspace(object):
    """
    Initializes a Molspace instance
    -------------------------------
    
    This class can be called via:
      * Full namespace syntax    : ``mooonpy.molspace.molspace.Molspace()``
      * Aliased namespace syntax : ``mooonpy.Molspace()``
    """

    def __init__(self, filename='', **kwargs):
        """        
        Initialization Parameters
        -------------------------
        filename : str, optional
            An optional filename to read and initialize a Molspace() instance
            with molecular system information (e.g. atoms, bonds, force field
            parameters ...). Supported file extensions:
                
              * LAMMPS datafile ``.data``
              * Tripos mol2 file ``.mol2``
              * SYBL mol file ``.mol``
        
            If no filename is provided the Molspace instance will be generated
            with no molecular system information.
            
            
        Attributes
        ----------
        N : int
            The order of the filter. For 'bandpass' and 'bandstop' filters,
            the resulting order of the final second-order sections ('sos')
            matrix is ``2*N``, with `N` the number of biquad sections
            of the desired system.
        Wn : array_like
            The critical frequency or frequencies. For lowpass and highpass
            filters, Wn is a scalar; for bandpass and bandstop filters,
            Wn is a length-2 sequence.
            
        Methods
        -------
        N : int
            The order of the filter. For 'bandpass' and 'bandstop' filters,
            the resulting order of the final second-order sections ('sos')
            matrix is ``2*N``, with `N` the number of biquad sections
            of the desired system.
        """

        # print(rcParams)

        # Get some basic config options from kwargs or setup defaults
        # print(kwargs)

        self.astyles = kwargs.pop('astyles', rcParams['molspace.astyles'])
        self.dsect = kwargs.pop('dsect', rcParams['molspace.read.dsect'])
        self.steps = kwargs.pop('steps', rcParams['molspace.read.steps'])

        # print(kwargs)

        # Build this object with some composition
        self.atoms: Atoms = Atoms(self.astyles)
        self.bonds: Bonds = Bonds()
        self.angles: Angles = Angles()
        self.dihedrals: Dihedrals = Dihedrals()
        self.impropers: Impropers = Impropers()
        self.clusters: Clusters = Clusters(self)
        self.ff: ForceField = ForceField()
        self.ptable: Ptable = Ptable()

        # Handle file initializations
        self.filename = filename
        self.header = ''
        if filename != '' and filename is not None:
            filename = Path(filename)
            if not bool(filename):  # check existence
                raise FileNotFoundError(f'{filename} was not found or is a directory')

            self.read_files(filename, dsect=self.dsect, steps=self.steps)
            if filename.endswith('.dump') and type(self.steps) != int:
                warnings.warn(
                    f'WARNING: Cannot init from a dump file with multiple steps specified. This returned the last step probably')

            # keys = self.bonds

    def copy(self):
        return deepcopy(self)
        ## This is absurdly slow, not sure why?

    def __str__(self):
        populated = []
        if self.bonds:
            populated.append('Bonds')
        if self.ff:
            populated.append('FF')
        ## Continue this for important attributes
        return f'Molspace with {len(self.atoms)} atoms and {populated} attributes)'

    def read_files(self, filename, dsect=['all'], steps=...):
        """
        Read files into Molspace
        Currently supports .data and .dump files.

        :param filename: path to file to read
        :type filename: str or Path
        :param dsect: Sections of data file to read, see mooonpy.molspace._files_io.read_lmp_data for more details
        :type dsect: list
        :param steps: Steps in dump file to read see mooonpy.molspace._files_io.read_lmp_dump for more details
        :type steps: list or int

        :return: Modifies self in-place if data or single step, returns dict if multiple steps
        :rtype: None or dict

        :Example:
            >>> from mooonpy import Molspace
            >>> MyPath = Path('Project/Monomers/DETDA.data')
            >>> mol = Molspace(MyPath, dsect=['all'])
            >>> MyPath2 = Path('Project/Monomers/DETDA.dump')
            >>> series = Molspace(MyPath2, steps=None)
        """
        root, ext = os.path.splitext(filename)
        if filename.endswith('.data'):
            if 'all' in dsect:
                dsect = ['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities']
            _files_io.read_lmp_data.read(self, filename, dsect)
            return None

        elif filename.endswith('.dump'):
            if steps is ...:  # none means all, ellipsis means to pull default from creation or rcParams
                steps = self.steps
            out_steps = _files_io.read_lmp_dump.read(filename, steps=steps, mol=self)
            if type(steps) is int: # if single, modify self with data
                self.atoms = out_steps.atoms
            return out_steps

        else:
            raise FileNotFoundError(f'{filename} was not found or is a directory')

    def write_files(self, filename, atom_style='full'):
        root, ext = os.path.splitext(filename)
        if filename.endswith('.data'):
            _files_io.write_lmp_data.write(self, filename, atom_style)
        if filename.endswith('.ff.script'):
            _files_io.write_lmp_ff_script.write(self, filename)

    def compute_pairs(self, cutoff, whitelist=None, blacklist=None, algorithm='DD_13', periodicity='ppp'):
        """
        Compute pairwise distances within cutoff between atoms
        """
        if algorithm == 'DD_13':
            domains, fractionals = domain_decomp_13(self.atoms, cutoff, whitelist, blacklist, periodicity)
            pairs = pairs_from_domains(self.atoms, cutoff, domains, fractionals)
        else:
            raise ValueError('Algorithm must be DD_13')

        return domains, pairs

    def compute_bond_length(self, periodicity='ppp'):
        return pairs_from_bonds(self.atoms, self.bonds, 'ppp')

    def compute_ADI(self):
        return ADI_from_bonds(self.bonds, self.angles, self.dihedrals, self.impropers)

    def compute_BADI_by_type(self, periodicity='ppp', comp_bond=True, comp_angle=True, comp_dihedral=True,
                             comp_improper=True):
        if comp_bond:
            pairs_from_bonds(self.atoms, self.bonds, periodicity)

        if comp_angle:
            angles = self.angles
        else:
            angles = None
        if comp_dihedral:
            dihedrals = self.dihedrals
        else:
            dihedrals = None
        if comp_improper:
            impropers = self.impropers
        else:
            impropers = None

        ADI_from_bonds(self.bonds, angles, dihedrals, impropers)
        bond_hist, angle_hist, dihedral_hist, improper_hist = BADI_by_type(self, self.ff.has_type_labels, comp_bond,
                                                                           comp_angle, comp_dihedral, comp_improper)
        return bond_hist, angle_hist, dihedral_hist, improper_hist

    def add_type_labels(self, labels):
        """
        Adds type labels for atom types and use 1st instance for each BADI type label

        could refactor the composition to a function somewhere
        """
        ff = self.ff
        ## Atom types
        for atom_type, coeff in ff.masses.items():
            if atom_type in labels:
                coeff.type_label = labels[atom_type]
                ff.pair_coeffs[atom_type].type_label = labels[atom_type]  # does this always exist
            else:
                raise Exception(f'Type label {atom_type} is not defined in input dictionary')

        ## Bond Types from example
        bond_types = {}
        for key, bond in self.bonds.items():
            if bond.type in bond_types: continue
            tl1 = labels[self.atoms[bond.ordered[0]].type]  # use load order for label
            tl2 = labels[self.atoms[bond.ordered[1]].type]
            bond_types[bond.type] = f"{tl1}-{tl2}"
            if len(bond_types) == len(ff.bond_coeffs): break  # found all examples
        # else: not every type has an example

        for type_, label in bond_types.items():
            ff.bond_coeffs[type_].type_label = label
            ## also cross terms?

        angle_types = {}
        for key, angle in self.angles.items():
            if angle.type in angle_types: continue
            tl1 = labels[self.atoms[angle.ordered[0]].type]
            tl2 = labels[self.atoms[angle.ordered[1]].type]
            tl3 = labels[self.atoms[angle.ordered[2]].type]
            angle_types[angle.type] = f"{tl1}-{tl2}-{tl3}"
            if len(angle_types) == len(ff.angle_coeffs): break

        for type_, label in angle_types.items():
            ff.angle_coeffs[type_].type_label = label

        dihedral_types = {}
        for key, dihedral in self.dihedrals.items():
            if dihedral.type in dihedral_types: continue
            tl1 = labels[self.atoms[dihedral.ordered[0]].type]
            tl2 = labels[self.atoms[dihedral.ordered[1]].type]
            tl3 = labels[self.atoms[dihedral.ordered[2]].type]
            tl4 = labels[self.atoms[dihedral.ordered[3]].type]
            dihedral_types[dihedral.type] = f"{tl1}-{tl2}-{tl3}-{tl4}"
            if len(dihedral_types) == len(ff.dihedral_coeffs): break

        for type_, label in dihedral_types.items():
            ff.dihedral_coeffs[type_].type_label = label

        improper_types = {}
        for key, improper in self.impropers.items():
            if improper.type in improper_types: continue
            tl1 = labels[self.atoms[improper.ordered[0]].type]
            tl2 = labels[self.atoms[improper.ordered[1]].type]
            tl3 = labels[self.atoms[improper.ordered[2]].type]
            tl4 = labels[self.atoms[improper.ordered[3]].type]
            improper_types[improper.type] = f"{tl1}-{tl2}-{tl3}-{tl4}"
            if len(improper_types) == len(ff.improper_coeffs): break

        for type_, label in improper_types.items():
            ff.improper_coeffs[type_].type_label = label

        ff.has_type_labels = True

    def generate_graph(self):
        return interface.generate_graph(self)

    def find_rings(self, ring_sizes: tuple[int] = (3, 4, 5, 6, 7)):
        """
        Finds all rings in the current Molspace instance. The size of rings searched for is
        set in the in the ring_sizes parameter.
        
        :Example:
        >>> import mooonpy
        >>> my_molecule = mooonpy.molspace('detda.mol2')
        >>> rings = my_molecule.find_rings(ring_sizes=(3,4,5,6,7))
        >>> print(rings)
        [(1, 2, 3, 4, 5, 6)]
        
        .. note::
            Each ring will be sorted in the order it was traversed in the graph and
            will be in the canonical form (e.g., canonical=(1,2,3) vs non_canonical=(2,3,1)).

        :param ring_sizes: Tuple containing ring sizes to search for in graph
        :type ring_sizes: tuple[int]
        :return: rings
        :rtype: list[tuple[int]]
        """
        graph = self.generate_graph()
        return ring_analysis.find_rings(graph, ring_sizes)

    def find_cumulative_neighs(self):
        return

    def update_elements(self, type2mass=None, type2element=None):
        """
        Updates every per-atom .element attribute with the element type for that atom.
        
        :Example:
        >>> import mooonpy
        >>> my_molecule = mooonpy.molspace('detda.data')
        >>> my_molecule.update_elements()
        
        .. note::
            This will also update the per-mass .element attribute in ff.masses dictionary as well.

        :param type2mass: Optional for setting the atom type to mass map (e.g. type2mass={1: 12.01115, 2: 1.008}). If not provided this will be generated from the masses section.
        :type type2mass: dict[int, float]
        :param type2element: Optional for setting the atom type to mass map (e.g. type2element={1: 'C', 2: 'H'}). If not provided this will be generated from the masses section and interally defined periodic table.
        :type type2element: dict[int, str]
        :return: None
        :rtype: None
        """
        if type2mass is None:
            type2mass = {i: self.ff.masses[i].coeffs[0] for i in self.ff.masses}
        if type2element is None:
            type2element = {i: self.ptable.mass2element(type2mass[i]) for i in type2mass}
        for atom_type, mass in self.ff.masses.items():
            mass.element = type2element[atom_type]
        for atom_id, atom in self.atoms.items():
            atom.element = type2element[atom.type]
        return

    def bonds_from_distances(self, periodicity: str = 'ppp'):
        """
        Resets the bonds in a molecular system, based in interatomic distances and valences of atoms. 
        The cutoff distances are set based on the summation of vdw radii per element involved in the
        bond. Each atom's valence is respected (e.g. the maximum number of allowable bonded atoms to
        a carbon atom is 4).
        
        :Example:
        >>> import mooonpy
        >>> my_molecule = mooonpy.molspace('detda.data')
        >>> my_molecule.bonds_from_distances(periodicity='fff')
        
        .. note::
            Before using this command the element per-atom info and element per-mass info need to be updated by
            running "my_molecule.update_elements()".

        :param periodicity: Optional for setting the periodicity of the model, where f=fixed and p=periodic. Each position is a face (e.g. periodicity='fpp' is fixed on X-faces and periodic in Y- and Z-faces)
        :type periodicity: str
        :return: None
        :rtype: None
        """
        find_bonds_from_distances(self, periodicity)
        return

    def remove_atoms(self, atom_ids):
        """ make smart and remove from some?"""
        atom_ids = set(atom_ids)
        for id_ in list(self.atoms.keys()):  # cant change size in dict loop, cache keys
            if id_ in atom_ids:  # set inclusion is O(1)
                del self.atoms[id_]
        for key in list(self.bonds.keys()):
            for key_id in key:
                if key_id in atom_ids:
                    del self.bonds[key]
                    break  # later ids would break
        for key in list(self.angles.keys()):
            for key_id in key:
                if key_id in atom_ids:
                    del self.angles[key]
                    break
        for key in list(self.dihedrals.keys()):
            for key_id in key:
                if key_id in atom_ids:
                    del self.dihedrals[key]
                    break
        for key in list(self.impropers.keys()):
            for key_id in key:
                if key_id in atom_ids:
                    del self.impropers[key]
                    break

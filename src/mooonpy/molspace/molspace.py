# -*- coding: utf-8 -*-
from . import _files_io as _files_io
from mooonpy.tools.file_utils import Path
from .atoms import Atoms
from .topology import Bonds, Angles, Dihedrals, Impropers
from .force_field import ForceField
from .distance import domain_decomp_13, pairs_from_bonds, pairs_from_domains, ADI_from_bonds, BADI_by_type
from .graph_theory import interface, ring_analysis
from mooonpy.rcsetup import rcParams
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

        # print(kwargs)

        # Build this object with some composition
        self.atoms: Atoms = Atoms(self.astyles)
        self.bonds: Bonds = Bonds()
        self.angles: Angles = Angles()
        self.dihedrals: Dihedrals = Dihedrals()
        self.impropers: Impropers = Impropers()
        self.ff: ForceField = ForceField()

        # Handle file initilaizations
        self.filename = filename
        self.header = ''
        if filename != '' and filename is not None:
            filename = Path(filename)
            if not bool(filename): # check existence
                raise FileNotFoundError(f'{filename} was not found or is a directory')

            self.read_files(filename, dsect=self.dsect)

            # keys = self.bonds

    def read_files(self, filename, dsect=['all']):
        root, ext = os.path.splitext(filename)
        if filename.endswith('.data'):
            if 'all' in dsect:
                dsect = ['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities']
            _files_io.read_lmp_data.read(self, filename, dsect)

        return None

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

    def compute_BADI_by_type(self,periodicity='ppp',comp_bond=True, comp_angle=True, comp_dihedral=True, comp_improper=True):
        if comp_bond:
            pairs_from_bonds(self.atoms, self.bonds, 'ppp')
        ADI_from_bonds(self.bonds, self.angles, self.dihedrals, self.impropers)
        bond_hist, angle_hist, dihedral_hist, improper_hist = BADI_by_type(self, self.ff.has_type_labels)
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
        
    def find_rings(self, ring_sizes: tuple[int]=(3,4,5,6,7)):
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

# -*- coding: utf-8 -*-
from mooonpy.tools.file_utils import smart_open
from mooonpy.tools.string_utils import string2digit


def read(mol, filename, sections):
    # Define sections to read (using inputs from user if they pass them)
    sections_mp: list[str] = ['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities', 'Ellipsoids', 'Lines',
                              'Triangles', 'Bodies']
    sections_tl: list[str] = ['Atom Type Labels', 'Bond Type Labels', 'Angle Type Labels', 'Dihedral Type Labels',
                              'Improper Type Labels']
    sections_ff: list[str] = ['Masses', 'Pair Coeffs', 'Bond Coeffs', 'Angle Coeffs', 'Dihedral Coeffs',
                              'Improper Coeffs', 'PairIJ Coeffs']
    sections_xt: list[str] = ['BondBond Coeffs', 'BondAngle Coeffs', 'AngleAngleTorsion Coeffs',
                              'EndBondTorsion Coeffs',
                              'MiddleBondTorsion Coeffs', 'BondBond13 Coeffs', 'AngleTorsion Coeffs',
                              'AngleAngle Coeffs']
    sections_fixes: list[str] = ['bond_react_props_internal']

    # Set up the sets that will be used for parsing
    sections_all: set[str] = set(sections_mp + sections_tl + sections_ff + sections_xt + sections_fixes)
    sections_coeffs: set[str] = sections_ff + sections_xt
    sections_kwargs: set[str] = set(sections)

    # Create shortcuts to all the generation factories for speed and ease of use
    atom_factory = mol.atoms.styles.atom_factory
    bond_factory = mol.bonds.bond_factory
    angle_factory = mol.angles.angle_factory
    dihedral_factory = mol.dihedrals.dihedral_factory
    improper_factory = mol.impropers.improper_factory
    coeffs_factory = mol.ff.coeffs_factory

    # Find "hard coded" atom style reads for performance and set atom_reader to slow and update later on
    hard_coded_atom_reading_styles = {i.split('_')[-1]: i for i in dir(mol.atoms.styles) if i.startswith('read_')}
    atom_reader = mol.atoms.styles.atom_fill

    # Open and read contents from file
    skip: int = 0
    section: str = ''
    ff_coeffs: None = None  # Will be a pointer to specifc ff_coeffs to update
    with smart_open(filename) as f:
        f = f.readlines()
        for n, string in enumerate(f):
            # skip line between section keywords and "top of the body"
            if skip > 0:
                skip -= 1
                continue

            # Toggle section "off" since a blank line will be at 
            # the "bottom of the body" (skip handles the blank line
            # at the "top of the body").
            elif string == '\n' or not string:
                section = ''
                continue

            # Deal with comments
            elif '#' in string:
                line = string.partition('#')
                data_str = line[0].strip()
                data_lst = data_str.split()
                comment = line[2].strip()
            else:
                data_str = string.strip()
                data_lst = data_str.split()
                comment = ''

            # -------------------------------------------------------#
            # Parse the computationally heavy parts 1st:
            #  - setting the section requires looking at each line
            #  - if we already know the section we can get extra
            #    performance by not having to set section flag for
            #    the entire stretch of the large sections
            # -------------------------------------------------------#
            if section == 'Atoms':
                atom = atom_factory()
                atom = atom_reader(atom, mol.atoms.style, data_lst)
                atom.comment = comment
                mol.atoms[atom.id] = atom

            elif section == 'Velocities':
                nid = int(data_lst[0])
                vx = float(data_lst[1])
                vy = float(data_lst[2])
                vz = float(data_lst[3])

                atom = mol.atoms[nid]
                atom.vx = vx
                atom.vy = vy
                atom.vz = vz

            elif section == 'Bonds':
                # nid = int(data_lst[0])
                type_id = string2digit(data_lst[1])  # This  could be a type label
                id1 = int(data_lst[2])
                id2 = int(data_lst[3])
                ordered = [id1, id2]

                bond = bond_factory()
                bond.ordered = ordered  # [id1, id2]
                bond.comment = comment
                bond.type = type_id
                key = tuple(sorted(ordered))
                mol.bonds[key] = bond

            elif section == 'Angles':
                # nid = int(data_lst[0])
                type_id = string2digit(data_lst[1])  # This  could be a type label
                id1 = int(data_lst[2])
                id2 = int(data_lst[3])
                id3 = int(data_lst[4])
                ordered = [id1, id2, id3]

                angle = angle_factory()
                angle.ordered = ordered  # [id1, id2, id3]
                angle.comment = comment
                angle.type = type_id
                if id1 < id3:
                    key = tuple(ordered)
                else:
                    key = (id3, id2, id1)
                mol.angles[key] = angle

            elif section == 'Dihedrals':
                # nid = int(data_lst[0])
                type_id = string2digit(data_lst[1])  # This  could be a type label
                id1 = int(data_lst[2])
                id2 = int(data_lst[3])
                id3 = int(data_lst[4])
                id4 = int(data_lst[5])
                ordered = [id1, id2, id3, id4]

                dihedral = dihedral_factory()
                dihedral.ordered = ordered  # [id1, id2, id3, id4]
                dihedral.comment = comment
                dihedral.type = type_id
                if id1 < id4:
                    key = tuple(ordered)
                else:
                    key = (id4, id3, id2, id1)
                mol.dihedrals[key] = dihedral

            elif section == 'Impropers':
                # nid = int(data_lst[0])
                type_id = string2digit(data_lst[1])  # This  could be a type label
                id1 = int(data_lst[2])
                id2 = int(data_lst[3])
                id3 = int(data_lst[4])
                id4 = int(data_lst[5])
                ordered = [id1, id2, id3, id4]

                improper = improper_factory()
                improper.ordered = ordered  # [id1, id2, id3, id4]
                improper.comment = comment
                improper.type = type_id
                outers = sorted([id1, id3, id4])
                key = (outers[0], id2, outers[1], outers[2])
                mol.impropers[key] = improper


            # Type labels can initialize a ff_coeffs build
            elif section in sections_tl and ff_coeffs is not None:
                typeID = int(data_lst[0])
                type_label = str(data_lst[1])

                params = coeffs_factory()
                params.comment = type_label
                params.type_label = type_label
                ff_coeffs[typeID] = params


            # Read-in force field parameters (type labels might have already initialized a 
            # ff_coeffs build - if not one will be initialized here)
            elif section in sections_coeffs and ff_coeffs is not None:
                # print(n, line, section)
                digits = [string2digit(string) for string in data_lst]
                typeID = digits[0]
                coeffs = digits[1:]
                
                # Some Coeffs section can be in the following format:
                #   Bond Coeffs  # morse class2
                #  
                #   1  morse     82.696   2.000     1.530          # c2          c3
                #   2  morse     82.696   2.000     1.501          # c2          cp
                #   3  class2    1.1010 345.000  -691.890  844.600 # c2        hpan
                #
                # In which case we will "strip" out the hydrid styles so read-in
                # coeffs are purely numeric and we will check for hybrid styles when
                # write the LAMMPS datafile or LAMMPS script file.
                if isinstance(coeffs[0], str):
                    style = coeffs[0]
                    del coeffs[0]
                else:
                    style = ff_coeffs.style
                                   
                # If typeID already exits, it means a type label built
                # this section, so only add what is needed in this case
                if typeID in ff_coeffs:
                    ff_coeffs[typeID].coeffs = coeffs
                    ff_coeffs[typeID].comment = comment
                    ff_coeffs[typeID].style = style
                    
                # Else a type label did not trigger a build of the coeffs
                # and build everything from scratch
                else:
                    # Build type label from comment, if type label
                    # doesnt already exist or set as read-in typeID
                    if comment:
                        type_label = '-'.join(comment.split())
                    else:
                        type_label = 'tl|{}'.format(str(typeID))

                    # Generate a params instance and add to it 
                    params = coeffs_factory(coeffs)
                    params.comment = comment
                    params.type_label = type_label
                    params.style = style
                    ff_coeffs[typeID] = params

                    # Get box dimensions
            elif 'xlo' in data_str and 'xhi' in data_str:
                mol.atoms.box.xlo = float(data_lst[0])
                mol.atoms.box.xhi = float(data_lst[1])
                continue
            elif 'ylo' in data_str and 'yhi' in data_str:
                mol.atoms.box.ylo = float(data_lst[0])
                mol.atoms.box.yhi = float(data_lst[1])
                continue
            elif 'zlo' in data_str and 'zhi' in data_str:
                mol.atoms.box.zlo = float(data_lst[0])
                mol.atoms.box.zhi = float(data_lst[1])
                continue
            elif 'xy' in data_str and 'xz' in data_str and 'yz' in data_str:
                mol.atoms.box.xy = float(data_lst[0])
                mol.atoms.box.xz = float(data_lst[1])
                mol.atoms.box.yz = float(data_lst[2])
                continue
            elif n == 0:
                mol.header = string
                continue


            # -----------------------------------------------------------------#
            # Toggle between sections. Toggling is expensive:                 
            #  - Requires each line to be checked, which means every line in
            #    Atoms, Bond, ... etc needs to be checked                     
            #  - If we use wise if/elif settings, once the Atoms, Bonds, ...  
            #    etc sections have been found we do not have to check         
            #    if data_str is a section to parse                            
            # -----------------------------------------------------------------#
            # elif data_str in sections_all:
            elif data_str[0].isalpha():
                skip = 1  # skip the line under each section keyword
                section = data_str
                if section not in sections_all:
                    raise Exception(f'ERROR {section} is not a supported LAMMPS datafile section')

                # Set flags for molecule data like atoms, bonds, etc ... Also check if user wants
                # that section read or not (if not set section to '', to skip reading that section)
                if section == 'Atoms':
                    mol.atoms.style = comment
                    if 'Atoms' not in sections_kwargs: section = ''

                    # Update atom reader to a quicker one (if supported)
                    if mol.atoms.style in hard_coded_atom_reading_styles:
                        reader_name = hard_coded_atom_reading_styles[mol.atoms.style]
                        atom_reader = getattr(mol.atoms.styles, reader_name)

                elif section == 'Bonds':
                    mol.bonds.style = comment
                    if 'Bonds' not in sections_kwargs:
                        section = ''
                elif section == 'Angles':
                    mol.angles.style = comment
                    if 'Angles' not in sections_kwargs:
                        section = ''
                elif section == 'Dihedrals':
                    mol.dihedrals.style = comment
                    if 'Dihedrals' not in sections_kwargs:
                        section = ''
                elif section == 'Impropers':
                    mol.impropers.style = comment
                    if 'Impropers' not in sections_kwargs:
                        section = ''
                elif section == 'Velocities' and 'Velocities' not in sections_kwargs:
                    section = ''

                # Type labels can initialize a ff dictionary (e.g. Atom Type
                # Labels, will generate the mol.ff.masses and then once masses
                # are ready, the coeffs will be updated at that point).
                elif section == 'Atom Type Labels':
                    ff_coeffs = mol.ff.masses
                    mol.ff.has_type_labels = True
                elif section == 'Bond Type Labels':
                    ff_coeffs = mol.ff.bond_coeffs
                    mol.ff.has_type_labels = True
                elif section == 'Angle Type Labels':
                    ff_coeffs = mol.ff.angle_coeffs
                    mol.ff.has_type_labels = True
                elif section == 'Dihedral Type Labels':
                    ff_coeffs = mol.ff.dihedral_coeffs
                    mol.ff.has_type_labels = True
                elif section == 'Improper Type Labels':
                    ff_coeffs = mol.ff.improper_coeffs
                    mol.ff.has_type_labels = True

                # Force field related parsing
                elif section == 'Masses':
                    ff_coeffs = mol.ff.masses
                    ff_coeffs.style = comment
                elif section == 'Pair Coeffs':
                    ff_coeffs = mol.ff.pair_coeffs
                    ff_coeffs.style = comment
                elif section == 'Bond Coeffs':
                    ff_coeffs = mol.ff.bond_coeffs
                    ff_coeffs.style = comment
                elif section == 'Angle Coeffs':
                    ff_coeffs = mol.ff.angle_coeffs
                    ff_coeffs.style = comment
                elif section == 'Dihedral Coeffs':
                    ff_coeffs = mol.ff.dihedral_coeffs
                    ff_coeffs.style = comment
                elif section == 'Improper Coeffs':
                    ff_coeffs = mol.ff.improper_coeffs
                    ff_coeffs.style = comment
                elif section == 'BondBond Coeffs':
                    ff_coeffs = mol.ff.bondbond_coeffs
                    ff_coeffs.style = comment
                elif section == 'BondAngle Coeffs':
                    ff_coeffs = mol.ff.bondangle_coeffs
                    ff_coeffs.style = comment
                elif section == 'AngleAngleTorsion Coeffs':
                    ff_coeffs = mol.ff.angleangletorsion_coeffs
                    ff_coeffs.style = comment
                elif section == 'EndBondTorsion Coeffs':
                    ff_coeffs = mol.ff.endbondtorsion_coeffs
                    ff_coeffs.style = comment
                elif section == 'MiddleBondTorsion Coeffs':
                    ff_coeffs = mol.ff.middlebondtorsion_coeffs
                    ff_coeffs.style = comment
                elif section == 'BondBond13 Coeffs':
                    ff_coeffs = mol.ff.bondbond13_coeffs
                    ff_coeffs.style = comment
                elif section == 'AngleTorsion Coeffs':
                    ff_coeffs = mol.ff.angletorsion_coeffs
                    ff_coeffs.style = comment
                elif section == 'AngleAngle Coeffs':
                    ff_coeffs = mol.ff.angleangle_coeffs
                    ff_coeffs.style = comment
                else:
                    # We need to toggle ff_coeffs between different sections
                    # so we do not add coeffs to the wrong dictionaries
                    ff_coeffs = None

    print_info = False
    # print_info = True
    if print_info:
        print('\n\n\nSTYLES:')
        print('header: ', mol.header[:50])
        print('atom style: ', mol.atoms.style)
        print('xbox: ', mol.atoms.box.xlo, mol.atoms.box.xhi)
        print('ybox: ', mol.atoms.box.ylo, mol.atoms.box.yhi)
        print('zbox: ', mol.atoms.box.zlo, mol.atoms.box.zhi)
        print('tilt: ', mol.atoms.box.xy, mol.atoms.box.xz, mol.atoms.box.yz)
        print('natoms: ', len(mol.atoms))
        print('nbonds: ', len(mol.bonds))
        print('nangles: ', len(mol.angles))
        print('ndihedrals: ', len(mol.dihedrals))
        print('nimpropers: ', len(mol.impropers))

        print('\n\nAtoms')
        for n, (key, value) in enumerate(mol.atoms.items()):
            if n < 5:
                print(key, value.type, value.x, value.y, value.z, value.comment, value.element, id(value))

        print('\n\nLINE GEN')
        line = mol.atoms.styles.atom_line(mol.atoms[1], style='full')
        print(line)

        print('\n\nBonds')
        for n, (key, value) in enumerate(mol.bonds.items()):
            if n < 5:
                print(key, value.type, value.ordered, value.bo, value.comment, id(value))

        print('\n\nAngles')
        for n, (key, value) in enumerate(mol.angles.items()):
            if n < 5:
                print(key, value.type, value.ordered, value.comment, id(value))

        print('\n\nDihedrals')
        for n, (key, value) in enumerate(mol.dihedrals.items()):
            if n < 5:
                print(key, value.type, value.ordered, value.comment, id(value))

        print('\n\nImpropers')
        for n, (key, value) in enumerate(mol.impropers.items()):
            if n < 5:
                print(key, value.type, value.ordered, value.comment, id(value))

        d = mol.ff.masses
        # d = mol.ff.pair_coeffs
        # d = mol.ff.bond_coeffs
        # d = mol.ff.angle_coeffs
        # d = mol.ff.dihedral_coeffs
        # d = mol.ff.improper_coeffs
        # d = mol.ff.bondbond_coeffs
        # d = mol.ff.bondangle_coeffs
        # d = mol.ff.angleangletorsion_coeffs
        # d = mol.ff.endbondtorsion_coeffs
        # d = mol.ff.middlebondtorsion_coeffs
        # d = mol.ff.bondbond13_coeffs
        # d = mol.ff.angletorsion_coeffs
        # d = mol.ff.angleangle_coeffs
        print('\n\nFF-CHECK: ', d.style)
        for key, value in d.items():
            print('{} {}   "{}"   "{}"'.format(key, value.coeffs, value.type_label, value.comment))

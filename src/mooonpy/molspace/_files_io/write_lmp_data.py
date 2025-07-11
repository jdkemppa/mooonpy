# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
June 26, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
# Function for stringing together float values for parameters
def string_parameters(coeff):
    string = ''
    for i in coeff:
        if isinstance(i, float):
            string += '{:>16.10f}'.format(i)
        elif isinstance(i, int):
            string += '{:>6}'.format(i)
        else:
            string += '{:>6}'.format(i)
    return string


# Function for writing lammps datafile
def write(mol, filename, atom_style):

    with open(filename,'w') as f:
        # Write header
        header = mol.header
        f.write(f'{header[-220:len(header)]}\n') # Make max header length of 220 characters

        # lines = ['line:{}\n'.format(i) for i in range(10)]
        # f.writelines(lines)

        # Write structure quantities
        if mol.atoms: f.write(f'{len(mol.atoms)} atoms\n')
        if mol.bonds: f.write(f'{len(mol.bonds)} bonds\n')
        if mol.angles: f.write(f'{len(mol.angles)} angles\n')
        if mol.dihedrals: f.write(f'{len(mol.dihedrals)} dihedrals\n')
        if mol.impropers: f.write(f'{len(mol.impropers)} impropers\n')
        f.write('\n')

        # Write structure quantity types
        if mol.ff.masses: f.write(f'{len(mol.ff.masses)} atom types\n')
        if mol.ff.bond_coeffs: f.write(f'{len(mol.ff.bond_coeffs)} bond types\n')
        elif mol.bonds:
            # added elif's 5-Jul-25 TDM. OVITO does not like no counts when bonds have types,
            # this is slow but will get a count if the ff is not populated
            bondtypes = set(bond.type for bond in mol.bonds.values())
            f.write(f'{len(bondtypes)} bond types\n')
        if mol.ff.angle_coeffs: f.write(f'{len(mol.ff.angle_coeffs)} angle types\n')
        elif mol.angles:
            angletypes = set(angle.type for angle in mol.angles.values())
            f.write(f'{len(angletypes)} angle types\n')
        if mol.ff.dihedral_coeffs: f.write(f'{len(mol.ff.dihedral_coeffs)} dihedral types\n')
        elif mol.dihedrals:
            dihedraltypes = set(dihedral.type for dihedral in mol.dihedrals.values())
            f.write(f'{len(dihedraltypes)} dihedral types\n')
        if mol.ff.improper_coeffs: f.write(f'{len(mol.ff.improper_coeffs)} improper types\n')
        elif mol.impropers:
            impropertypes = set(improper.type for improper in mol.impropers.values())
            f.write(f'{len(impropertypes)} improper types\n')
        f.write('\n')

        # write box size
        f.write('{:^16.10f} {:^16.10f} {:^5} {:5}\n'.format(mol.atoms.box.xlo, mol.atoms.box.xhi, 'xlo', 'xhi'))
        f.write('{:^16.10f} {:^16.10f} {:^5} {:5}\n'.format(mol.atoms.box.ylo, mol.atoms.box.yhi, 'ylo', 'yhi'))
        f.write('{:^16.10f} {:^16.10f} {:^5} {:5}\n'.format(mol.atoms.box.zlo, mol.atoms.box.zhi, 'zlo', 'zhi'))
        if any([mol.atoms.box.xy, mol.atoms.box.xz, mol.atoms.box.yz]):
            f.write('{:^16.10f} {:^16.10f} {:^16.10f} {} {} {}\n'.format(mol.atoms.box.xy, mol.atoms.box.xz, mol.atoms.box.yz, 'xy', 'xz', 'yz'))

        # Write Atom Type Labels if user wants
        if mol.ff.has_type_labels and mol.ff.masses:
            f.write('\nAtom Type Labels\n\n')
            type_ids = sorted(mol.ff.masses.keys())
            used = set()
            for i in type_ids:
                coeffs = mol.ff.masses[i]
                type_label = coeffs.type_label
                f.write('{:^3} {:^2}\n'.format(i, type_label))
                if ' ' in type_label:
                    print('WARNING Atom Type Label: {} has been white space.'.format(type_label))
                    print(__file__)
                if type_label in used:
                    print('WARNING Atom Type Label: {} has been defined multiple times.'.format(type_label))
                    print(__file__)
                used.add(type_label)

        # Write Bond Type Labels if user wants
        if mol.ff.has_type_labels and mol.ff.bond_coeffs:
            f.write('\nBond Type Labels\n\n')
            type_ids = sorted(mol.ff.bond_coeffs.keys())
            used = set()
            for i in type_ids:
                coeffs = mol.ff.bond_coeffs[i]
                type_label = coeffs.type_label
                f.write('{:^3} {:^2}\n'.format(i, type_label))
                if ' ' in type_label:
                    print('WARNING Bond Type Label: {} has been white space.'.format(type_label))
                    print(__file__)
                if type_label in used:
                    print('WARNING Bond Type Label: {} has been defined multiple times.'.format(type_label))
                    print(__file__)
                used.add(type_label)

        # Write Angle Type Labels if user wants
        if mol.ff.has_type_labels and mol.ff.angle_coeffs:
            f.write('\nAngle Type Labels\n\n')
            type_ids = sorted(mol.ff.angle_coeffs.keys())
            used = set()
            for i in type_ids:
                coeffs = mol.ff.angle_coeffs[i]
                type_label = coeffs.type_label
                f.write('{:^3} {:^2}\n'.format(i, type_label))
                if ' ' in type_label:
                    print('WARNING Angle Type Label: {} has been white space.'.format(type_label))
                    print(__file__)
                if type_label in used:
                    print('WARNING Angle Type Label: {} has been defined multiple times.'.format(type_label))
                    print(__file__)
                used.add(type_label)

        # Write Dihedral Type Labels if user wants
        if mol.ff.has_type_labels and mol.ff.dihedral_coeffs:
            f.write('\nDihedral Type Labels\n\n')
            type_ids = sorted(mol.ff.dihedral_coeffs.keys())
            used = set()
            for i in type_ids:
                coeffs = mol.ff.dihedral_coeffs[i]
                type_label = coeffs.type_label
                f.write('{:^3} {:^2}\n'.format(i, type_label))
                if ' ' in type_label:
                    print('WARNING Dihedral Type Label: {} has been white space.'.format(type_label))
                    print(__file__)
                if type_label in used:
                    print('WARNING Dihedral Type Label: {} has been defined multiple times.'.format(type_label))
                    print(__file__)
                used.add(type_label)

        # Write Improper Type Labels if user wants
        if mol.ff.has_type_labels and mol.ff.improper_coeffs:
            f.write('\nImproper Type Labels\n\n')
            type_ids = sorted(mol.ff.improper_coeffs.keys())
            used = set()
            for i in type_ids:
                coeffs = mol.ff.improper_coeffs[i]
                type_label = coeffs.type_label
                f.write('{:^3} {:^2}\n'.format(i, type_label))
                if ' ' in type_label:
                    print('WARNING Improper Type Label: {} has been white space.'.format(type_label))
                    print(__file__)
                if type_label in used:
                    print('WARNING Improper Type Label: {} has been defined multiple times.'.format(type_label))
                    print(__file__)
                used.add(type_label)

        # Write massses
        f.write('\nMasses\n\n')
        type_ids = sorted(mol.ff.masses.keys())
        for i in type_ids:
            coeff = mol.ff.masses[i]
            parms = coeff.coeffs
            if coeff.comment:
                comment = '# {}'.format(coeff.comment)
            else:
                comment = ''
            f.write('{:^3} {} {}\n'.format(i, string_parameters(parms), comment))

        # Write all sections that have the "Coeffs ending"
        potentials = [i for i in dir(mol.ff) if i.endswith('coeffs')]
        write_order = ['pair_coeffs', 'bond_coeffs', 'angle_coeffs', 'dihedral_coeffs', 'improper_coeffs',
                       'bondbond_coeffs', 'bondangle_coeffs', 'angleangletorsion_coeffs', 'endbondtorsion_coeffs',
                       'middlebondtorsion_coeffs', 'bondbond13_coeffs', 'angletorsion_coeffs', 'angleangle_coeffs']
        for attr in write_order:
            if attr not in potentials:
                print('WARNING potential: {} is not being written to {}.'.format(attr, filename))
                print(__file__)
                continue
            
            potential = getattr(mol.ff, attr)
            if not potential: continue
            
            if potential.style:
                style_hint = '# {}'.format(potential.style)
            else: style_hint = ''
            
            type_ids = sorted(potential.keys())
            if not type_ids: continue
            f.write('\n{} {}\n\n'.format(potential.keyword, style_hint))
            line_styles = mol.ff.get_per_line_styles(attr)
            unique_line_styles = set(line_styles.values())
            line_lengths = set()
            for i in type_ids:
                coeff = potential[i]
                
                if len(unique_line_styles) == 1:
                    style = ''
                else: style = '{:^8}'.format(coeff.style)
                
                if coeff.comment:
                    comment = '# {}'.format(coeff.comment)
                else: comment = ''
                
                line = '{:^3} {} {}'.format(i, style, string_parameters(coeff.coeffs))
                line_lengths.add(len(line))                
                f.write('{:<{s}} {}\n'.format(line, comment, s=max(line_lengths)))

        # Write atoms and velocities
        if mol.atoms:
            if mol.atoms.style:
                style_hint = '# {}'.format(mol.atoms.style)
            else: style_hint = ''
            f.write('\nAtoms # {}\n\n'.format(atom_style))

            # Setup the atom line generation script
            hard_coded_atom_line_styles = {i.split('_')[-1]:i for i in dir(mol.atoms.styles) if i.startswith('line_')}
            atom_line = mol.atoms.styles.atom_line
            if atom_style in hard_coded_atom_line_styles:
                line_name = hard_coded_atom_line_styles[atom_style]
                atom_line = getattr(mol.atoms.styles, line_name)

            # Finally write the atoms
            atoms = sorted(mol.atoms.keys())
            for i in atoms:
                atom = mol.atoms[i]

                if atom.comment:
                    comment = '# {}'.format(atom.comment)
                else: comment = ''

                line = atom_line(atom, style=atom_style)
                f.write('{} {}\n'.format(line, comment))

            f.write('\nVelocities\n\n')
            for i in atoms:
                atom = mol.atoms[i]
                f.write('{:^6} {:^16.10f} {:^16.10f} {:^16.10f}\n'.format(i, atom.vx, atom.vy, atom.vz))

        # Write bonds
        if mol.bonds:
            f.write('\nBonds\n\n')
            for i, keys in enumerate(mol.bonds, 1):
                topo = mol.bonds[keys]

                if topo.comment:
                    comment = '# {}'.format(topo.comment)
                else: comment = ''

                id1, id2 = topo.ordered
                f.write('{:^6} {:^3} {:^6} {:^6} {}\n'.format(i, topo.type, id1, id2, comment))

        # Write angles
        if mol.angles:
            f.write('\nAngles\n\n')
            for i, keys in enumerate(mol.angles, 1):
                topo = mol.angles[keys]

                if topo.comment:
                    comment = '# {}'.format(topo.comment)
                else: comment = ''

                id1, id2, id3 = topo.ordered
                f.write('{:^6} {:^3} {:^6} {:^6} {:^6} {}\n'.format(i, topo.type, id1, id2, id3, comment))

        # Write dihedrals
        if mol.dihedrals:
            f.write('\nDihedrals\n\n')
            for i, keys in enumerate(mol.dihedrals, 1):
                topo = mol.dihedrals[keys]

                if topo.comment:
                    comment = '# {}'.format(topo.comment)
                else: comment = ''

                id1, id2, id3, id4 = topo.ordered
                f.write('{:^6} {:^3} {:^6} {:^6} {:^6} {:^6} {}\n'.format(i, topo.type, id1, id2, id3, id4, comment))

        # Write impropers
        if mol.impropers:
            f.write('\nImpropers\n\n')
            for i, keys in enumerate(mol.impropers, 1):
                topo = mol.impropers[keys]

                if topo.comment:
                    comment = '# {}'.format(topo.comment)
                else: comment = ''

                id1, id2, id3, id4 = topo.ordered
                f.write('{:^6} {:^3} {:^6} {:^6} {:^6} {:^6} {}\n'.format(i, topo.type, id1, id2, id3, id4, comment))

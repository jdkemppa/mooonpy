# -*- coding: utf-8 -*-
import warnings
from typing import Union, Optional
from mooonpy.tools.file_utils import smart_open, Path
from mooonpy.tools.string_utils import string2digit

_Molspace = None

def _get_molspace_class():
    """
    Lazy import on first run
    """
    global _Molspace
    if _Molspace is None:
        from mooonpy.molspace.molspace import Molspace
        _Molspace = Molspace
    return _Molspace

def read(filename: Union[str, Path], steps: Optional[int|list]=None, mol=None):
    """
    Read dump file into Molspace.
    steps is a whitelist of steps to read. An int will return the first matching molspace,
    an iterable (set, list, dict) will return all matching Molspaces in a {step:mol} dict,
    None will return all Molspaces in a {step:mol} dict.

    If mol is not none, each step will modify a copy of the input molspace. None will start with an empty molspace.
    The atoms attribute is cleared from the copy, so all other attributes are copied

    :param filename: Path to dump file
    :type filename: str or Path
    :param steps: Whitelist of steps to read
    :type steps: int or iterable of int
    :param mol: Base Molspace
    :type mol: Molspace or None
    """
    if mol is None:
        Molspace = _get_molspace_class() # circular so cannot use in top import
        # from mooonpy.molspace.molspace import Molspace # slow for many calls
        mol = Molspace() # empty
    # if mol is None:
    #     from mooonpy.molspace.molspace import Molspace
    #     mol = Molspace()
    atom_factory = mol.atoms.styles.atom_factory
    style = mol.atoms.styles

    section: str = ''
    box_info: str = ''
    box_line: int = 0
    fields: list = []
    step = None
    if type(steps) is int:
        current_mol = None
        out_mols = None
    else:
        current_mol = None
        out_mols = {}

    with smart_open(filename) as f:
        # f = f.readlines()
        # for n, string in enumerate(f):
        for string in f:
            ## Mode switching
            if string.startswith('ITEM'):
                if 'TIMESTEP' in string:
                    if steps is None:
                        pass
                    elif isinstance(steps, int) and section == 'ATOMS':
                        # end of single step read
                        break
                    elif not isinstance(steps, int) and len(out_mols) == len(steps) and section == 'ATOMS':
                        break  # all found
                    section = 'TIMESTEP'
                elif section == 'SKIP':  # disables below elifs for moving past quick
                    continue
                elif 'NUMBER' in string:
                    section = 'NUMBER'
                elif 'BOX' in string:
                    section = 'BOX'
                    box_info = string.strip().split()[3:]  # not currently stored anywhere
                    box_line = 0
                elif 'ATOMS' in string:
                    section = 'ATOMS'
                    fields = string.strip().split()[2:]
                    functions = [style.func_aliases[field] for field in fields]
                    # print(fields,functions)
                else:
                    raise Exception('Unknown section {}'.format(string))

                continue  # do not run sections

            elif section == 'ATOMS':  # goes first because called many times
                atom_info = string.strip().split()
                try:
                    atom_dict = {field: fun(info) for field, fun, info in zip(fields, functions, atom_info)}
                except:
                    print(fields, functions, atom_info)
                    raise Exception('Unknown atom info {}'.format(atom_info))
                id_ = atom_dict['id']  # not removed and is added to atom object
                # id_ = atom_dict.pop('id')

                if id_ in mol.atoms:
                    atom = mol.atoms[id_].copy()  # inherit base info not in dump
                else:
                    atom = atom_factory()

                for field, value in atom_dict.items():
                    try:
                        setattr(atom, field, value)
                    except:
                        raise Exception('Unknown field {}. Add to Styles defaults'.format(field))
                current_mol.atoms[id_] = atom

            elif section == 'SKIP':
                continue

            elif section == 'TIMESTEP':
                step = int(string)
                # print(step, steps)
                if steps is None:  # all mode
                    current_mol = mol.copy()
                    out_mols[step] = current_mol
                elif isinstance(steps, int) and steps == step:  # single mode True
                    current_mol = mol.copy()
                    current_mol.atoms.clear()
                    out_mols = current_mol  # pointer
                elif isinstance(steps, int) and steps != step:  # single mode False
                    section = 'SKIP'
                    current_mol = None
                elif step not in steps:  # whitelist False
                    section = 'SKIP'
                elif step in steps:  # whitelist True
                    current_mol = mol.copy()
                    current_mol.atoms.clear()
                    out_mols[step] = current_mol
                else:
                    raise Exception('Control Flow error in string: {}'.format(string))
                continue

            elif section == 'BOX':
                box = current_mol.atoms.box  # get pointer
                # box_info contains periodicity, do we need an attribute for that?
                box_line += 1  # x,y,z
                line = string.strip().split()
                if box_line == 1:
                    box.xlo = string2digit(line[0])
                    box.xhi = string2digit(line[1])
                    if len(line) == 3:
                        box.xy = string2digit(line[2])  # triclinics (untested)
                elif box_line == 2:
                    box.ylo = string2digit(line[0])
                    box.yhi = string2digit(line[1])
                    if len(line) == 3:
                        box.xz = string2digit(line[2])
                elif box_line == 3:
                    box.zlo = string2digit(line[0])
                    box.zhi = string2digit(line[1])
                    if len(line) == 3:
                        box.yz = string2digit(line[2])
                else:
                    raise Exception('Unknown box line {}'.format(string))

            elif section == 'NUMBER':
                continue  # nothing to do with this info?
            else:
                raise Exception('Unknown section {}'.format(section))

        if step is None:
            warnings.warn('No Steps Match Selection')
            return None
    return out_mols


if __name__ == '__main__':
    from mooonpy.molspace.molspace import Molspace
    from mooonpy.molspace.atom_styles import Styles
    from mooonpy.molspace._files_io.read_lmp_dump import read
    ## Local path to example files
    test_dir = Path(__file__).dir() / '..\\..\\..\\..\\examples\\EPON_862\\lmp_dump'

    ## create extra fields
    ## This is best used for fixes. the box fraction shown here is now a default, but was tested before
    # Styles.styles['box'] = ('xu', 'yu', 'zu')  # adding styles before first Molspace call is safe
    # Styles.defaults.update({'xu': 0.0, 'yu': 0.0, 'zu': 0.0})
    # Styles.func_aliases.update({'xu': float, 'yu': float, 'zu': float})

    ## get PCFF model
    bonded = Molspace(test_dir / 'small_epon.data')
    print(bonded.atoms[1].type == 10) # true thing to check
    start_x = bonded.atoms[1].x
    ## Call function directly
    dump_dict = read(test_dir / 'small_epon_densify.dump', steps=[0, 1000], mol=bonded)
    ## ^ mol can inherit BADI or None makes empty

    ## __init__ from specified timestep returns Molspace, not dict, same as data file
    one_step = Molspace(test_dir / 'small_epon_densify.dump',steps=1000)
    ## ^ mol=___ inheritance does not work here. Could fix by piping args through init. but there are other ways to solve the problem

    ## This does not work correctly, because "init"-ing from multiple files does not know what to return
    bad_many_step = Molspace(test_dir / 'small_epon_densify.dump')
    ## ^ AGAIN, DO NOT USE THIS
    ## ^ I think it's using the last step in single mode
    ## ^ raises a warning now

    ## This does work to init one or many molspaces as a dict
    all = Molspace().read_files(test_dir / 'small_epon_densify.dump')
    ## ^ behavior should be the same the "dump_dict" example with steps the same

    ## Dict of steps inheriting BADI from bonded
    double = bonded.read_files(test_dir / 'small_epon_densify.dump', steps=[0, 1000])
    ## single inherit from bonded
    single = bonded.read_files(test_dir / 'small_epon_densify.dump', steps=1000)
    ## ^ note this does NOT modify bonded
    ## ^ It's worth noting that the deepcopy used in each step is very slow for some reason

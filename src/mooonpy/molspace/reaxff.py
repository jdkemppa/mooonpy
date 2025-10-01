# -*- coding: utf-8 -*-
from collections import defaultdict
from typing import Optional
from mooonpy.tools.file_utils import Path, smart_open
from math import sqrt, exp


class ReaxFF:
    def __init__(self, topofile=None):
        self.name: str = topofile
        self.comment: str = ''
        self.genreral = {}
        self.atom_param_list = []  # the FF file does not actually use these and are missing several parameters
        self.bond_param_list = []  # so these lines are actually useless
        self.element_param = {}
        self.bond_param = {}
        self.ele_list = []
        if topofile is not None:
            self.read(topofile)

    def read(self, topofile: Optional[Path | str]):
        topofile = Path(topofile)
        self.name = str(topofile.basename())

        section_rows = 0
        row_counter = 0
        # param_index = 0
        element = None

        with smart_open(topofile) as f:
            mode = 'comment'
            for ii, string in enumerate(f):
                string = string.strip()
                if '!' in string:
                    if 'general parameters' in string:
                        mode = 'general'
                        continue
                    elif 'Nr of atoms' in string:
                        mode = 'atoms'
                        self.atom_param_list += string.partition(';')[2].split(';')[1:]
                        # ignore ist delimiter, and atomID section
                        section_rows = 1
                        continue

                    elif 'Nr of bonds' in string:
                        mode = 'bonds'
                        self.bond_param_list += string.partition(';')[2].split(';')[2:]
                        # ignore ist delimiter, and at1 at2
                        section_rows = 1
                        continue

                    elif mode == 'general':  # has ! as delimiter
                        splits = string.partition('!')
                        self.genreral[splits[2].strip()] = float(splits[0].strip())

                    else:
                        mode = 'none'

                elif mode == 'atoms':
                    if ';' in string:
                        self.atom_param_list += string.split(';')
                        section_rows += 1
                    else:
                        splits = string.split()
                        row_counter += 1
                        if row_counter == 1:
                            element = splits.pop(0)
                            # print(element)
                            self.element_param[element] = {}
                            self.ele_list.append(element)
                            self.element_param[element]['lgre'] = 0.0  # defaults
                            self.element_param[element]['lgcij'] = 0.0
                            # param_index = 0
                            keys = ['r_s', 'valency', 'mass', 'r_vdw', 'epsilon', 'gamma', 'r_pi', 'valency_e']
                        elif row_counter == 2:
                            keys = ['alpha', 'gamma_w', 'valency_boc', 'p_ovun5', 'n.u.', 'chi', 'eta', 'p_hbond']
                        elif row_counter == 3:
                            keys = ['r_pi_pi', 'p_lp2', 'n.u.', 'b_o_131', 'b_o_132', 'b_o_133', 'bcut_acks2', 'n.u.']
                        elif row_counter == 4:
                            keys = ['p_ovun2', 'p_val3', 'n.u.', 'valency_val', 'p_val5', 'rcore2', 'ecore2', 'acore2']
                        elif row_counter == 5:  # not always there?
                            keys = ['lgcij', 'lgre']
                        else:
                            continue
                        # print(len(keys))
                        for key, value in zip(keys, splits):
                            if key == 'n.u.': continue
                            # print(element, key, value)
                            self.element_param[element][key] = float(value)

                        if row_counter == section_rows:
                            row_counter = 0  # reset for new element next line


                elif mode == 'bonds':
                    if ';' in string:
                        self.bond_param_list += string.split(';')
                        section_rows += 1
                    else:
                        splits = string.split()
                        row_counter += 1
                        if row_counter == 1:
                            element_1 = self.ele_list[int(splits.pop(0)) - 1]
                            element_2 = self.ele_list[int(splits.pop(0)) - 1]
                            param1 = self.element_param[element_1]
                            param2 = self.element_param[element_2]
                            # print(element_1, element_2)
                            bond_key = (element_1, element_2)
                            self.bond_param[bond_key] = {}
                            ## make composite attributes. reaxff_ffield.cpp lines 359-380
                            self.bond_param[bond_key]['r_s'] = 0.5 * (param1['r_s'] + param2['r_s'])
                            if param1['r_pi'] > 0 and param2['r_pi'] > 0:
                                self.bond_param[bond_key]['r_p'] = 0.5 * (param1['r_pi'] + param2['r_pi'])
                            else:
                                self.bond_param[bond_key]['r_p'] = None  # diable for H, F and other single bonds
                            if param1['r_pi_pi'] > 0 and param2['r_pi_pi'] > 0:
                                self.bond_param[bond_key]['r_pp'] = 0.5 * (param1['r_pi_pi'] + param2['r_pi_pi'])
                            else:
                                self.bond_param[bond_key]['r_pp'] = None
                            self.bond_param[bond_key]['p_boc3'] = sqrt(param1['b_o_132'] * param2['b_o_132'])
                            self.bond_param[bond_key]['p_boc4'] = sqrt(param1['b_o_131'] * param2['b_o_131'])
                            self.bond_param[bond_key]['p_boc5'] = sqrt(param1['b_o_133'] * param2['b_o_133'])
                            self.bond_param[bond_key]['D'] = sqrt(param1['epsilon'] * param2['epsilon'])
                            self.bond_param[bond_key]['alpha'] = sqrt(param1['alpha'] * param2['alpha'])
                            self.bond_param[bond_key]['r_vdW'] = 2.0 * sqrt(param1['r_vdw'] * param2['r_vdw'])
                            self.bond_param[bond_key]['gamma_w'] = sqrt(param1['gamma_w'] * param2['gamma_w'])
                            self.bond_param[bond_key]['gamma'] = pow(param1['gamma'] * param2['gamma'], -1.5)

                            self.bond_param[bond_key]['rcore'] = sqrt(param1['rcore2'] * param2['rcore2'])
                            self.bond_param[bond_key]['ecore'] = sqrt(param1['ecore2'] * param2['ecore2'])
                            self.bond_param[bond_key]['acore'] = sqrt(param1['acore2'] * param2['acore2'])

                            keys = ['De_s', 'De_p', 'De_pp', 'p_be1', 'p_bo5', 'v13cor', 'p_bo6', 'p_ovun1']
                        elif row_counter == 2:
                            keys = ['p_be2', 'p_bo3', 'p_bo4', 'n.u.', 'p_bo1', 'p_bo2', 'ovc', 'n.u.']
                            # no 8th, is comment wrong?

                        for key, value in zip(keys, splits):
                            if key == 'n.u.': continue
                            # print(bond_key, key, value)
                            self.bond_param[bond_key][key] = float(value)

                        if row_counter == section_rows:
                            # print(self.bond_param[bond_key]['r_s'], self.bond_param[bond_key]['r_p'],
                            #       self.bond_param[bond_key]['r_pp'])
                            # print(self.bond_param[bond_key]['De_s'], self.bond_param[bond_key]['De_p'],
                            #       self.bond_param[bond_key]['De_pp'])
                            row_counter = 0  # reset for new element next line


                elif mode == 'comment':
                    self.comment = string

    def sort_ele(self, in_list):
        """
        Usage:
        bond_key = tuple(reax.sort_ele([atom1.element, atom2.element]))
        """
        indexes = [self.ele_list.index(ele) for ele in in_list]
        indexes.sort()
        out_list = [self.ele_list[ii] for ii in indexes]
        return out_list

    def pair_bo_prime(self, ele1, ele2, dist):
        bond_key = self.sort_ele([ele1, ele2])
        params = self.bond_param[tuple(bond_key)]

        BO_s_prime = exp(params['p_bo1'] * pow((dist / params['r_s']), params['p_bo2']))
        if params['r_p'] is not None:
            BO_p_prime = exp(params['p_bo3'] * pow((dist / params['r_p']), params['p_bo4']))
        else:
            BO_p_prime = 0.0
        if params['r_pp'] is not None:
            BO_pp_prime = exp(params['p_bo5'] * pow((dist / params['r_pp']), params['p_bo6']))
        else:
            BO_pp_prime = 0.0

        return BO_s_prime, BO_p_prime, BO_pp_prime

    def BO_average(self, series, base=None, remove_single=False, remove_H2=False, steps=None, cutoff=3, periodicity='ppp' ):
        """

        """
        import mooonpy as mp
        if base is None:
            base = mp.Molspace()
        elif isinstance(base, mp.Molspace):
            pass
        elif isinstance(base, str) or isinstance(base, Path):
            base = mp.Molspace(base)
        else:
            raise TypeError('base is not a Molspace or Path')
        bond_factory = base.bonds.bond_factory
        base.update_elements()


        if isinstance(series, mp.Molspace):
            series = {0:series}
        elif isinstance(series, dict):
            pass
        elif isinstance(series, str) or isinstance(series, Path):
            series = mp.Molspace().read_files(series, steps=steps)
        else:
            raise TypeError('series is not a Molspace or Path')

        new_bonds = defaultdict(list)
        for ss, dump in series.items():
            dump.atoms.wrap()
            domains, pairs = dump.compute_pairs(cutoff, periodicity=periodicity)
            for key, pair in pairs.items():
                atom1 = dump.atoms[key[0]]
                atom2 = dump.atoms[key[1]]
                ## This will be corrected BO eventually
                BO_prime = self.pair_bo_prime(atom1.element, atom2.element, pair.distance)
                new_bonds[key].append(sum(BO_prime))

        out_mol = series[list(series.keys())[0]]
        for id_, atom in out_mol.atoms.items():
            if id_ in base.atoms:
                atom.type = base.atoms[id_].type

        bond_avg = {}
        for key, BO_list in new_bonds.items():
            if len(BO_list) < len(series): continue  # ignore pass through cutoff
            BO_avg = sum(BO_list) / len(BO_list)
            if BO_avg < 0.3: continue
            ## could also use min in the range
            bond_avg[key] = BO_avg

        bond_avg = dict(sorted(bond_avg.items(), key=lambda item: item[1], reverse=True))
        ## ^ high to low BO
        atom_graph = {id_: [] for id_ in out_mol.atoms.keys()}
        atom_total = {id_: 0.0 for id_ in out_mol.atoms.keys()} # BO total per atom

        for key, BO in bond_avg.items():
            atom1 = out_mol.atoms[key[0]]
            atom2 = out_mol.atoms[key[1]]
            val1 = self.element_param[atom1.element]['valency']
            val2 = self.element_param[atom2.element]['valency']

            if remove_H2 and atom1.element == 'H' and atom2.element == 'H': continue  # removed if it cannot form other bond
            # if atom_total[key[0]] > val1-BO or atom_total[key[1]] > val2-BO: continue # over valenced after this bond

            if len(atom_graph[key[0]]) < val1 and len(atom_graph[key[1]]) < val2: # not over bonded
                atom_graph[key[0]].append(key[1])
                atom_graph[key[1]].append(key[0])
                atom_total[key[0]] += BO # use attribute in atom object instead?
                atom_total[key[1]] += BO
                out_mol.bonds[key] = bond_factory()
                out_mol.bonds[key].type = 1  # dummy
                out_mol.bonds[key].ordered = key  # needed for write. change so fallback to key?

        if remove_single:
            for id_, graph in atom_graph.items():
                if len(graph) == 0:
                    out_mol.atoms.pop(id_)  # remove lone atoms

        return out_mol
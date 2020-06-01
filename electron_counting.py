# -*- coding: utf-8 -*

from __future__ import absolute_import, division, print_function, unicode_literals
import multiprocessing
import os
import numpy as np
import matplotlib.pyplot as plt
from progress.bar import Bar
from scipy.optimize import minimize
from itertools import combinations
from ase.visualize import view
from ase.db import connect
import random
from ase.build import sort
from ase import Atom
from ase.io import read, write

seed = 1234
random.seed(seed)
np.random.seed(seed)

db = connect('./qm9.db')

ucfc_cutoff = {'H': 1.65, 'C':1.86, 'O':1.86, 'F':1.85, 'N':1.83}
bond_cutoff = {'H': 1.40, 'C':2.20, 'O':2.20, 'F':2.20, 'N':2.0}

vlnc_cnt = {'H': 1, 'C':4, 'O':6, 'F':7, 'N':5}
full_cnt = {'H': 2, 'C':8, 'O':8, 'F':8, 'N':8}

def tag_ucfc(atoms, del_idx):
    ret = []
    for i in range(len(atoms)):
        if i == del_idx:
            ret.append(-1) # del
        elif atoms.get_distance(i, del_idx) < ucfc_cutoff[atoms[del_idx].symbol]:
            ret.append(0) # UC
        else:
            ret.append(1) # FC
    return ret

def check_ucfc_by_electron_counting(atoms, del_idx, ucfc_tag_lst):
    electron_cnt = [0] * len(atoms)
    for i in range(len(atoms)):
        electron_cnt[i] += vlnc_cnt[atoms[i].symbol]
    while(True):
        modified = False
        for i in range(len(atoms)):
            for j in range(i+1, len(atoms)):
                if i==del_idx or j==del_idx:
                    continue
                # elif atoms.get_distance(i, j) < (bond_cutoff[atoms[i].symbol] + bond_cutoff[atoms[j].symbol])/2 and \
                elif atoms.get_distance(i, j) < (ucfc_cutoff[atoms[i].symbol] + ucfc_cutoff[atoms[j].symbol])/2 and \
                        electron_cnt[i] < full_cnt[atoms[i].symbol] and \
                        electron_cnt[j] < full_cnt[atoms[j].symbol]:
                    electron_cnt[i] += 1
                    electron_cnt[j] += 1
                    modified = True
                else:
                    pass
        if modified == False:
            break
    err_cnt = 0
    for i in range(len(atoms)):
        if i == del_idx:
            continue
        if electron_cnt[i] != full_cnt[atoms[i].symbol] and ucfc_tag_lst[i] == 1:
            err_cnt += 1
        elif electron_cnt[i] == full_cnt[atoms[i].symbol] and ucfc_tag_lst[i] == 0:
            err_cnt += 1
        else:
            pass
    # if err_cnt != 0:
    #     print('-' * 100)
    #     print(electron_cnt, ucfc_tag_lst, err_cnt)
    #     print(del_idx)
    #     view(atoms)
    return err_cnt

if __name__ == '__main__':
    err_cnt = 0
    total_cnt = 0
    bar = Bar('checking', max=db.count())
    for row in db.select():
        bar.next()
        atoms = row.toatoms()
        del_atom_idx = random.randint(0, len(atoms)-1)
        tag_lst = tag_ucfc(atoms, del_atom_idx)
        total_cnt += len(tag_lst)
        err_cnt += check_ucfc_by_electron_counting(atoms, del_atom_idx, tag_lst)
        # if err_cnt != 0:
        #     break
        # if row.id % 100 == 0:
        #     print('error count, total count: ', err_cnt, total_cnt)
    bar.finish()
    print('error count, total count: ', err_cnt, total_cnt)
    pass

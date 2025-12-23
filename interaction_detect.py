# -*- coding: utf-8 -*-
# @Time    : 2025/6/18 14:10
# @Author  : xiang.zhang
# @FileName: temp.py
# @Software: PyCharm
# @usage   : identify interaction between two chains with cif structure file, generate pkl file, which need to be further processed with interaction_detect.py
import pickle
import sys

import pandas as pd
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 180)
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pymol import cmd
from tqdm import tqdm

import multiprocessing as mp


def init_pymol(model_path):
    cmd.reinitialize()  # 重置PyMOL环境
    # cmd.fetch('1qox')  # 示例：加载PDB 1QOX
    cmd.load(model_path)
    cmd.remove('solvent')  # 移除水分子
    cmd.hide('everything')
    cmd.show('cartoon')  # 显示卡通结构
    cmd.util.cbc()  # 按链着色

def find_contact_custom(selection, selection2, mode=2, sele_name='t1', sele_name2='t2', cutoff=3.6, idx=1):
    init_pymol(model_path)
    # print(f"{selection} and {selection2} in")
    all_don_acc1 = selection + " and (donor or acceptor)"
    all_don_acc2 = selection2 + " and  (donor or acceptor)"
    contact_result = {}
    s1 = cmd.select(sele_name, all_don_acc1)
    s2 = cmd.select(sele_name2, all_don_acc2)
    if not (s1 and s2):
        return contact_result
    pbar = tqdm(cmd.get_model(sele_name).atom, ncols=150, desc=f'{selection}::{selection2}::{idx}', ascii=False, position=idx, leave=False)
    for atom_s in pbar:
        cmd.select('s', f'index {atom_s.index}')
        for atom_q in cmd.get_model(sele_name2).atom:
            cmd.select('q_t', f'index {atom_q.index}')
            distance = cmd.distance('d', 's', 'q_t', mode=mode, cutoff=cutoff)
            distance_big = cmd.distance('d_b', 's', 'q_t', mode=mode, cutoff=cutoff)
            cmd.delete('q_t')
            cmd.delete('d')
            cmd.delete('d_b')
            if distance == 0 or distance >= cutoff:
                continue
            if f"{atom_s.index}-{atom_q.index}" in contact_result.keys() or f"{atom_q.index}-{atom_s.index}" in contact_result.keys():
                continue
            contact_result[f"{atom_s.index}-{atom_q.index}"] = distance
        cmd.delete('s')
        pbar.update(1)
    pbar.close()
    cmd.delete(sele_name)
    cmd.delete(sele_name2)
    contact_result['param'] = {
        "selection": selection,
        "selection2": selection2,
        "mode":mode,
        "sele_name":sele_name,
        "sele_name2":sele_name2,
        "cutoff":cutoff
    }
    return contact_result

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', default=None, type=str)
    parser.add_argument('--cutoff', default=3.6, type=float)
    parser.add_argument('--mode', default=2, type=int, choices=[0,2])
    parser.add_argument('--selection', default='Chain A', type=str)
    parser.add_argument('--selection2', default='chain B', type=str)
    parser.add_argument('--prefix', default='Cas', type=str)
    args = parser.parse_args()
    mode_desc = {
        0: "all",
        2: "polar"
    }
    selection = args.selection
    # selection = 'Chain A'
    selection2 = args.selection2
    # selection2 = 'Chain B'
    mode = args.mode
    # mode=2
    cutoff = args.cutoff
    # cutoff=3.6
    prefix = args.prefix

    model_path = args.pdb
    if not model_path.endswith('cif'):
        sys.exit('cif format pdb file required')
    cif_dict = MMCIF2Dict(model_path)
    cif_info = {
        "model_num": cif_dict['_atom_site.pdbx_PDB_model_num'],
        "atom_level": cif_dict['_atom_site.group_PDB'],
        "atom_idx": cif_dict['_atom_site.id'],
        "seq_id": cif_dict['_atom_site.label_seq_id'],
        "seq_id_vis": cif_dict['_atom_site.auth_seq_id'],
        "chain_id": cif_dict['_atom_site.label_asym_id'],
        "atom_symbol": cif_dict['_atom_site.type_symbol'],
        "atom_label": cif_dict['_atom_site.label_atom_id'],
        "residue_name": cif_dict['_atom_site.label_comp_id']
    }

    df = pd.DataFrame(cif_info)

    print(f'parameters: prefix: {prefix}\nselection: {selection}\nselection2: {selection2}\nmode: {mode}\ncutoff: {cutoff}')
    l = df[df.chain_id == 'A'].seq_id.astype(int).max()
    parm_list = list(map(lambda x: [str(x), str(x + 9)], list(range(1, l + 1, 10))))
    parm_list[-1][-1] = str(l)
    for idx, i in enumerate(parm_list):
        parm_list[idx] = [f"{selection} and resi {i[0]}-{i[1]}", f"{selection2}"]
    for idx, i in enumerate(parm_list):
        for p in [mode, 't1', 't2', cutoff]:
            i += [p]
        i += [idx + 1]
    assert mode in mode_desc.keys()
    processor = mp.Pool(10, initializer=tqdm.set_lock, initargs=(mp.RLock(),))
    result = processor.starmap(find_contact_custom, parm_list)
    with open(f'{prefix}_result_{selection2}_{mode_desc[mode]}_cutoff{cutoff}.pkl','wb') as h:
        pickle.dump(result, h)

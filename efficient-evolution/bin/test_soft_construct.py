# -*- coding: utf-8 -*-
# @Time    : 2024/4/7 13:38
# @Author  : xiang.zhang
# @FileName: test_soft_construct.py
# @Software: PyCharm
# @usage   :
import os, sys
import pandas as pd
from pathlib import Path
#from amis import reconstruct_multi_models
sys.path.extend(['/home/zhangxiang/gitFolder/efficient-evolution'])
from bin.amis import reconstruct_multi_models
from bin.amis import get_model_name, diff, reconstruct, soft_reconstruct

gene = 'mAb114-VH'
sequence = 'EVQLVESGGGLIQPGGSLRLSCAASGFALRMYDMHWVRQTIDKRLEWVSAVGPSGDTYYADSVKGRFAVSRENAKNSLSLQMNSLTAGDTAIYYCVRSDRGVAGLFDSWGQGILVTVSS'
alpha = None

model_names  = ['/data/.data2/esm_data/esm1b_t33_650M_UR50S.pt',
 '/data/.data2/esm_data/esm1v_t33_650M_UR90S_1.pt',
 '/data/.data2/esm_data/esm1v_t33_650M_UR90S_2.pt',
 '/data/.data2/esm_data/esm1v_t33_650M_UR90S_3.pt',
 '/data/.data2/esm_data/esm1v_t33_650M_UR90S_4.pt',
 '/data/.data2/esm_data/esm1v_t33_650M_UR90S_5.pt']

mutations_models_origin, mutations_model_names_origin = reconstruct_multi_models(
    sequence,
    model_names,
    alpha=alpha,
    return_names=True
)


mutations_models_soft, mutations_model_names_soft = reconstruct_multi_models(
    sequence,
    model_names,
    alpha=1,
    return_names=True
)

mutations_models_soft.keys().__len__()
mutations_models_origin.keys().__len__()

common_mutations = list(set(mutations_models_soft.keys()).intersection(mutations_models_origin.keys()))

for cm in common_mutations:
    print(mutations_models_origin.get(cm) == mutations_models_soft.get(cm))



mutation_info = pd.DataFrame()
for mut, info in sorted(mutations_model_names_origin.items(), key=lambda item: item[0]):
    pos, wt_aa, mut_aa = mut
    for model_info in info:
        model_name, (wt_prob, mut_prob) = model_info
        rec = pd.DataFrame({
            "geneIdentifier":[gene],
            "model":[Path(model_name).stem],
            "pos": [pos],
            "wt_aa": [wt_aa],
            "mut_aa": [mut_aa],
            "wt_prob": [wt_prob],
            "mut_prob": [mut_prob],
            "mut_seq":[sequence[:pos ]+mut_aa+sequence[pos + 1: ]],
            "wt_mut_compare":[f"{sequence[:pos ]}({wt_aa}/{mut_aa}){sequence[pos + 1: ]}"]
        })
        mutation_info = mutation_info._append(rec)

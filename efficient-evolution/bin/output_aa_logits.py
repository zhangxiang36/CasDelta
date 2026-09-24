# -*- coding: utf-8 -*-
# @Time    : 2025/10/16 9:14
# @Author  : xiang.zhang
# @FileName: output_aa_logits.py
# @Software: PyCharm
# @usage   :
import os, sys
import pandas as pd
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 180)
import numpy as np
from pathlib import Path
#from amis import reconstruct_multi_models
sys.path.extend(['/home/zhangxiang/gitFolder/efficient-evolution'])
from bin.amis import reconstruct_multi_models
import scipy
from bin.amis import get_model_name

alpha=1.0
# Cas delta-1 HM
seq_name = "CasDelta1-HM"
sequence = "MSKDGDKMKHSKQEVSEDKGVKKDPSVFRCQNLAILEMREEVADYYADLQADYRHYLPLMILWGQNGFMEPDETGLNFSAVLPSTAKTRCENRYGAQIAADLKLIADDLPSSLYAIFTNFYQPLTGVKRVGDRWKFVTSLNTVSDEEFAKKVKPLLRDYKGVTEEQWRRFREERVAGKTADQVWDDLDSILRNPTLAALKKIHDKHMFRVKRPESKLYCDALAHMLSSSYDSWKKINEARQLELQAMRDEQQGLLANKTTNSLYNLVCEFADELESQNYGLTRRIFFIARKEWEDEKSHYLPVIKLLVTEKYLPLLQADNKTVSPAVKAWDICSTLKQRSPYIRFVPPTKDYPIPFGQTGRGHKFTLKEEAGKVCFLVPGLTKDEMPLSGSHYFSSLRIDKVGTDFKISFRHRRKIRRKKSARVTVKQPRINGTVTEVGILRRNGKFFLRVAYKIAYPQHNIDLKSYFASAAPSTQKLELLPAQMRVAGVDLNISDPIVVSKADIFKGEDGGPLSVLDYGSGKVVGEPVIICKDTSRAGKISSLATKCRTLRDVIRNFRHSLHQGTPLPSKSLEFLQEVPVAQEDLKSPSPRYLIQTWIKYIRTEMKRLHYKVWLDGYCHTSEAMRMLDLVDQVAFLMKSYENIHVKPGHRIPKKTTASLKIRETSRTRFRRHVSRLVGARVVRACEDCQMIFLEDLSTGFREGSNNSLRRLFSAGQLKNSIVMAAEKVGKAVVFVCKDGTSIKDPVFGLHGLRPSKKTPGIPPDKKRLYVRRDGKIGYINADQAASLNVLLVGLSHSVAKYKFFVKDGKLQDDEEVAAQQGKKPRKIGRRLERLLKMNFDTPARLFFSLEKDQVVPCYEATETPYTGWVYVYHDRLLTTRQRDQQVNQLAQEVKDLLRDGVKIPKWDVIPDCDNSCLGFSPCLVLDDETAAVSVV"
model_names = [os.path.join('/data/.data2/esm_data/', i) for i in
               [
                    'esm1b_t33_650M_UR50S.pt',
                    'esm1v_t33_650M_UR90S_1.pt',
                    'esm1v_t33_650M_UR90S_2.pt',
                    'esm1v_t33_650M_UR90S_3.pt',
                    'esm1v_t33_650M_UR90S_4.pt',
                    'esm1v_t33_650M_UR90S_5.pt'
                ]
               ]
# mutations_models, mutations_model_names = reconstruct_multi_models(
#     sequence,
#     model_names,
#     alpha,
#     return_names=True
# )

# model.alphabet_.tok_to_idx[seq[pos]]
exclude = {'B', 'J', 'O', 'U', 'X', 'Z', '-', '.'}
writer_probs = pd.ExcelWriter(f"{seq_name}_probs.xlsx")
writer_probsCompare = pd.ExcelWriter(f"{seq_name}_probsCompare.xlsx")

for model_n in model_names:
    model = get_model_name(model_n)
    model_name = model_n.split('/')[-1].replace('.pt','')
    print(model_name)
    logits = model.predict_sequence_prob(sequence)
    probs = scipy.special.softmax(logits, axis=1)
    aa_to_idx = model.alphabet_.tok_to_idx
    idx_to_aa = {item:key for key, item in aa_to_idx.items()}
    df_probs = pd.DataFrame(probs)
    df_probs.index = ["<cls>"]+list(sequence) + ["<eos>"]
    df_probs.columns = aa_to_idx.keys()
    un_select = list(exclude) + model.alphabet_.all_special_tokens + ['<null_1>']
    df_probs = df_probs.loc[~df_probs.index.isin(un_select), ~df_probs.columns.isin(un_select)].T
    df_probs.to_excel(writer_probs, sheet_name=model_name)
    for col_idx in range(df_probs.shape[1]):
        wt_aa = df_probs.columns[col_idx]
        wt_prob = df_probs.iloc[:,col_idx].loc[wt_aa]
        df_probs.iloc[:, col_idx] = np.log(df_probs.iloc[:, col_idx]/wt_prob)
    df_probs.to_excel(writer_probsCompare, sheet_name=model_name)

writer_probs.close()
writer_probsCompare.close()

writer_probs = pd.ExcelWriter(f"{seq_name}_probs_masked.xlsx")
writer_probsCompare = pd.ExcelWriter(f"{seq_name}_probsCompare_masked.xlsx")

from tqdm import tqdm
for model_n in model_names:
    model = get_model_name(model_n)
    model_name = model_n.split('/')[-1].replace('.pt','')
    print(model_name)
    outputs = []
    for idx, aa in enumerate(tqdm(sequence)):
        s_idx = idx+1
        sequence_mask = sequence[:idx] + '<mask>' + sequence[(idx+1):]
        logits = model.predict_sequence_prob(sequence_mask)
        probs = scipy.special.softmax(logits, axis=1)
        outputs.append(probs[s_idx,:])
    comp = np.zeros(shape=(33,), dtype=np.float32)
    outputs = [comp] + outputs + [comp]
    outputs[0][model.alphabet_.tok_to_idx['<cls>']] = 1.
    outputs[-1][model.alphabet_.tok_to_idx['<eos>']] = 1.
    probs = np.vstack(outputs)
    aa_to_idx = model.alphabet_.tok_to_idx
    idx_to_aa = {item: key for key, item in aa_to_idx.items()}
    df_probs = pd.DataFrame(probs)
    df_probs.index = ["<cls>"] + list(sequence) + ["<eos>"]
    df_probs.columns = aa_to_idx.keys()
    un_select = list(exclude) + model.alphabet_.all_special_tokens + ['<null_1>']
    df_probs = df_probs.loc[~df_probs.index.isin(un_select), ~df_probs.columns.isin(un_select)].T
    df_probs.to_excel(writer_probs, sheet_name=model_name)
    for col_idx in range(df_probs.shape[1]):
        wt_aa = df_probs.columns[col_idx]
        wt_prob = df_probs.iloc[:, col_idx].loc[wt_aa]
        df_probs.iloc[:, col_idx] = np.log(df_probs.iloc[:, col_idx] / wt_prob)
    df_probs.to_excel(writer_probsCompare, sheet_name=model_name)

writer_probs.close()
writer_probsCompare.close()



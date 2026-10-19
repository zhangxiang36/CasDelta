# -*- coding: utf-8 -*-
# @Time    : 2026/10/14 15:48
# @Author  : xiang.zhang
# @FileName: tnp.py
# @Software: PyCharm
# @usage   :
import numpy as np
import pandas as pd
import re

def allele_desc_sum(allele_desc):
    result = {}
    for p in allele_desc.split(';'):
        p_del, p_is_wt, p_perecent = re.findall('(^.+?)-(True|False)-(.+?$)', p)[0]
        if not result.get(int(p_del), 0):
            result[int(p_del)] = 0
        result[int(p_del)] += float(p_perecent)
    return result

df = pd.read_excel(r"SangerICE_NGScompare.xlsx")

results = pd.DataFrame()
for rec in df[['Event','allele_type','allele_type_NGS','allele_desc','allele_desc_NGS']].itertuples(index=False, name=None):
    ev,allele_type,allele_type_NGS,allele_desc,allele_desc_NGS = rec
    ev_allele_desc_ICE = {key:value*100 for key, value in allele_desc_sum(allele_desc).items()}
    ev_allele_desc_NGS = allele_desc_sum(allele_desc_NGS)
    d1 = pd.DataFrame([ev_allele_desc_ICE]).T
    d1['sample_id'] = ev
    d1['del_len'] = d1.index.tolist()
    d1.rename(columns={0:'AlleleFrequency'}, inplace=True)
    d2 = pd.DataFrame([ev_allele_desc_NGS]).T
    d2['sample_id'] = ev
    d2['del_len'] = d2.index.tolist()
    d2.rename(columns={0:'AlleleFrequency'}, inplace=True)
    d = pd.merge(d1, d2, on=['sample_id', 'del_len'], how='outer', suffixes=('_ICE', '_NGS'))
    d = d.iloc[:, [1,2,0,3]]
    results = pd.concat([results, d])

results_common = results.dropna(axis=0, how='any')

results_overall = pd.DataFrame()
for rec in df[['Event','allele_type','allele_type_NGS','allele_desc','allele_desc_NGS']].itertuples(index=False, name=None):
    ev, allele_type, allele_type_NGS, allele_desc, allele_desc_NGS = rec
    ev_allele_desc_ICE = {key: value * 100 for key, value in allele_desc_sum(allele_desc).items()}
    ev_allele_desc_NGS = allele_desc_sum(allele_desc_NGS)
    indel_pct_ICE = sum([value for key, value in ev_allele_desc_ICE.items() if key != 0])
    indel_pct_NGS = sum([value for key, value in ev_allele_desc_NGS.items() if key != 0])
    r = pd.DataFrame(
        {
            'sample_id': [ev],
            'indel_pct_ICE': [indel_pct_ICE],
            'indel_pct_NGS': [indel_pct_NGS],
            "allele_type_ICE": [allele_type],
            "allele_type_NGS": [allele_type_NGS],
        }
    )
    results_overall = pd.concat([results_overall, r])

print("Pearson correlation coefficient between Indel percents of sanger-ICE and NGS is: ", np.corrcoef(results_overall['indel_pct_ICE'], results_overall['indel_pct_NGS'])[0][1])
print(f'NGS find {results.dropna(subset=['AlleleFrequency_NGS']).shape[0]} edited alleles, in which {results_common.shape[0]} edited alleles are also found with Sanger-ICE pipeline')
print(f'Among 12 samples tested, {(results_overall.allele_type_ICE == results_overall.allele_type_NGS).sum()} samples have the same allele type (het, chimeric, etc.)')

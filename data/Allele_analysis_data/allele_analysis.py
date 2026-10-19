# -*- coding: utf-8 -*-
# @Time    : 2026/10/13 15:43
# @Author  : xiang.zhang
# @FileName: allele_analysis.py
# @Software: PyCharm
# @usage   : Allele frequency extraction from ICE output (*.contribs.json)

import pandas as pd
import copy, re
import argparse
import json
def judgeAllele(content, filter_threshold = 0.05, wt_threshold = 0.05):
    """

    :param content: dict parsed from *.contribs.json
    :return:
    """
    content_copy = copy.deepcopy(content)
    alleles = []
    contribs_list = [i for i in content_copy['contribs_list']]
    for idx, contrib in enumerate(contribs_list):
        if contrib['rel_abundance'] < filter_threshold:
            continue
        contrib['indel'] = contrib['indel']['total']
        if contrib['wt'] and contrib['rel_abundance'] > wt_threshold:
            if not alleles:
                alleles = [contrib]
            else:
                alleles.insert(0,contrib)
        else:
            alleles.append(contrib)
    df_sample = pd.DataFrame(alleles)

    wt = sum(df_sample.wt)
    if wt:
        if df_sample.indel.abs().sum() > 0:
            allele = 'mutant'
            allele_num = df_sample[~df_sample.wt].indel.unique().__len__()
            allele_type = 'het' if allele_num == 1 else 'chimeric'
        else:
            allele = 'WT'
            allele_type = 'WT'
    else:
        allele = 'mutant'
        allele_num = df_sample[~df_sample.wt].indel.unique().__len__()
        if allele_num == 1:
            allele_type = 'homo'
        elif allele_num == 2:
            allele_type = 'bi-allele'
        elif allele_num > 2:
            allele_type = 'chimeric'
        else:
            allele_type = 'unknown'
    indel0_num = df_sample[df_sample.indel == 0].shape[0]
    if indel0_num > 0:
        if any(~df_sample[df_sample.indel==0].wt):
            has_snp = True
        else:
            has_snp = False
    else:
        has_snp = False

    assert( allele is not None )
    assert( allele_type is not None)
    allele_str = []
    for rec in df_sample[['indel', 'wt', 'rel_abundance']].to_records(index=False):
        allele_str.append('-'.join(map(lambda x: str(x), rec)))

    return (allele, allele_type, has_snp, df_sample.rel_abundance.sum().round(2), ';'.join(allele_str))


def allele_stat(allele_desc):
    """

    :param allele_desc: str, e.g.  '-28-False-0.16;-27-False-0.16;-1-False-0.09'
    :return: e.g. [28, 27, 1]
    """
    if "0-True" in allele_desc:
        allele_desc = re.sub(';{0,1}0-True-\\d\\.\\d{1,2};{0,1}','', allele_desc)
    mutation_type = []
    for mut in allele_desc.split(';'):
        units = mut.split('-')
        if mut.startswith('-'):
            units[0] = '-'
            al = ''.join(units[:2])
        else:
            al = units[0]
        mutation_type.append(
            int(al)
        )
    return list(set(mutation_type))

def main():
    parser = argparse.ArgumentParser(description='Allele frequency extraction from ICE output (*.contribs.json)')
    parser.add_argument('--input', '-i', help='ICE output file (*.contribs.json) to be processed', required=True)
    parser.add_argument('--filter_threshold',
                        help='Allele with percentage lower than threshold will be excluded, default 0.05 (5 percent)',
                        default=0.05,
                        type=float)
    args = parser.parse_args()

    input_file = args.input
    filter_threshold = args.filter_threshold

    with open(input_file, 'r') as handle:
        content = json.load(handle)

    allele_note, allele_type, has_snp, all_abundance, allele_desc = judgeAllele(content, filter_threshold)
    mutation_type = allele_stat(allele_desc)
    print(allele_note, allele_type, allele_desc, ';'.join(map(str, mutation_type)))

if __name__ == '__main__':
    main()

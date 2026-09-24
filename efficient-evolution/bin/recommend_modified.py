import os, sys
import pandas as pd
from pathlib import Path
# this script was from recommend.py, with some modifications for screening candidates with alpha and K parameters conviniently.
sys.path.extend(['/home/zhangxiang/gitFolder/CasDelta/efficient-evolution']) # modify if needed
from bin.amis import reconstruct_multi_models
# usage:
# python bin/recommend_modified.py --gene Cas12j_19 --alpha 1 --sequence PSYKSSRVLVRDVPEELVDHYERSHRVAAFFMRLLLAMRREPYSLRMRDGTEREVDLDETDDFLRSAGCEEPDAVSDDLRSFALAVLHQDNPKKRAFLESENCVSILCLEKSASGTRYYKRPGYQLLKKAIEEEWGWDKFEASLLDERTGEVAEKFAALSMEDWRRFFAARDPDDLGRELLKTDTREGMAAALRLRERGVFPVSVPEHLDLDSLKAAMASAAERLKSWLACNQRAVDEKSELRKRFEEALDGVDPEKYALFEKFAAELQQADYNVTKKLVLAVSAKFPATEPSEFKRGVEILKEDGYKPLWEDFRELGFVYLAERKWERRRGGAAVTLCDADDSPIKVRFGLTGRGRKFVLSAAGSRFLITVKLPCGDVGLTAVPSRYFWNPSVGRTTSNSFRIEFTKRTTENRRYVGEVKEIGLVRQRGRYYFFIDYNFDPEEVSDETKVGRAFFRAPLNESRPKPKDKLTVMGIDLGINPAFAFAVCTLGECQDGIRSPVAKMEDVSFDSTGLRGGIGSQKLHREMHNLSDRCFYGARYIRLSKKLRDRGALNDIEARLLEEKYIPGFRIVHIEDADERRRTVGRTVKEIKQEYKRIRHQFYLRYHTSKRDRTELISAEYFRMLFLVKNLRNLLKSWNRYHWTTGDRERRGGNPDELKSYVRYYNNLRMDTLKKLTCAIVRTAKEHGATLVAMENIQRVDRDDEVKRRKENSLLSLWAPGMVLERVEQELKNEGILAWEVDPRHTSQTSCITDEFGYRSLVAKDTFYFEQDRKIHRIDADVNAAINIARRFLTRYRSLTQLWASLLDDGRYLVNVTRQHERAYLELQTGAPAATLNPTAEASYELVGLSPEEEELAQTRIKRKKREPFYRHEGVWLTREKHREQVHELRNQVLALGNAKIPEIRT

def parse_args():
    esmv1_model_names = [os.path.join('/data/.data2/esm_data/', i) for i in \
                    [
                        'esm1b_t33_650M_UR50S.pt',
                        'esm1v_t33_650M_UR90S_1.pt',
                        'esm1v_t33_650M_UR90S_2.pt',
                        'esm1v_t33_650M_UR90S_3.pt',
                        'esm1v_t33_650M_UR90S_4.pt',
                        'esm1v_t33_650M_UR90S_5.pt'
                    ]
                ]
    esmv2_model_names = [os.path.join('/data/.data2/esm_data/checkpoints/', i) for i in \
                    [
                         # "esm2_t6_8M_UR50D.pt",
                         # "esm2_t12_35M_UR50D.pt",
                         # "esm2_t30_150M_UR50D.pt",
                         "esm2_t33_650M_UR50D.pt",
                         "esm2_t36_3B_UR50D.pt",
                         "esm2_t48_15B_UR50D.pt"
                    ]
                ]
    import argparse
    parser = argparse.ArgumentParser(
        description='Recommend substitutions to a wildtype sequence'
    )
    parser.add_argument('-s','--sequence', type=str,
                        help='Wildtype sequence')
    parser.add_argument(
        '--model-names',
        type=str,
        default=  # [ 'esm1b', 'esm1v1', 'esm1v2', 'esm1v3', 'esm1v4', 'esm1v5', ]
            # [os.path.join('/backup/.data2/esm_data/checkpoints/',i) for i in \
            [os.path.join('/data/.data2/esm_data/checkpoints/', i) for i in \
                [
                     "esm2_t6_8M_UR50D.pt",
                     "esm2_t12_35M_UR50D.pt",
                     "esm2_t30_150M_UR50D.pt",
                     "esm2_t33_650M_UR50D.pt",
                     "esm2_t36_3B_UR50D.pt",
                     "esm2_t48_15B_UR50D.pt"
                ]
            ],
        nargs='+',
        help='Type of language model (e.g., esm1b, esm1v1)'
    )
    parser.add_argument(
        '--model-version',
        type=str,
        default=  # [ 'esm1b', 'esm1v1', 'esm1v2', 'esm1v3', 'esm1v4', 'esm1v5', ]
        # [os.path.join('/backup/.data2/esm_data/checkpoints/',i) for i in \
        "esmv2",
        choices=['esmv1','esmv2','all'],
        required=True,
        help='Version of ESM language model (e.g., esmv1, esmv2, or both of them)'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=1.,
        help='alpha stringency parameter'
    )
    parser.add_argument('--topK', type=int, default=2, help='K stringency parameter, see paper for details')
    parser.add_argument(
        '--cuda',
        type=str,
        default='cuda',
        help='cuda device to use'
    )
    parser.add_argument(
        '--gene',
        type=str,
        default=None,
        help='gene name or identifier'
    )
    args = parser.parse_args()
    args.models = {'esmv1':esmv1_model_names, "esmv2":esmv2_model_names, 'all': esmv1_model_names+esmv2_model_names}
    return args

if __name__ == '__main__':
    args = parse_args()
    K = args.topK
    if ":" in args.cuda:
        os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda.split(':')[-1]

    mutations_models, mutations_model_names = reconstruct_multi_models(
        args.sequence,
        args.models[args.model_version],
        alpha=args.alpha,
        return_names=True
    )

    if not args.alpha:
        import sys
        sys.exit()

    mutation_info = pd.DataFrame()
    gene = 'Not provided' if not args.gene else args.gene
    print(f"Gene: {gene} is processing")
    for mut, info in sorted(mutations_model_names.items(), key=lambda item: item[0]):
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
                "mut_seq":[args.sequence[:pos ]+mut_aa+args.sequence[pos + 1: ]],
                "wt_mut_compare":[f"{args.sequence[:pos ]}({wt_aa}/{mut_aa}){args.sequence[pos + 1: ]}"]
            })
            mutation_info = mutation_info._append(rec)
    mutation_info.sort_values(['model','pos','mut_prob'], inplace=True, ascending=[True, True, False])
    mutation_info.reset_index(drop=True, inplace=True)
    output_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '{}_{}_ef-evo.xlsx'.format(gene, args.model_version)
    )
    mutation_info['improve'] = mutation_info.mut_prob / mutation_info.wt_prob
    # df.groupby(['pos', 'wt_aa'], as_index=False)['improve'].rank('max', ascending=False)
    mutation_info['rank'] = mutation_info.groupby(['pos', 'wt_aa'], as_index=False)['mut_prob'].rank('first', ascending=False)
    # df['max_exists'] = df.groupby(['pos','wt_aa','mut_aa'])['rank'].transform(lambda x:(x == 1).any())
    mutation_info['max_exists'] = mutation_info['rank'] == 1
    mutation_info[mutation_info.max_exists].drop_duplicates(subset=['pos', 'wt_aa', 'mut_aa'])
    mutation_info['model_num'] = mutation_info.groupby(['pos', 'wt_aa', 'mut_aa'])['model'].transform(lambda x: len(x))
    df_candidates = mutation_info[mutation_info.max_exists & (mutation_info.model_num >= K)].drop_duplicates(
        subset=['pos', 'wt_aa', 'mut_aa']).sort_values('pos')
    df_candidates.drop(['rank', 'max_exists'], axis='columns', inplace=True)
    df_candidates = df_candidates.iloc[:, [0, 1, 2, 3, 4, 5, 6, 9, 10, 7, 8]]
    writer = pd.ExcelWriter(output_path)
    df_candidates.to_excel(writer, index=False, sheet_name='candidates')
    mutation_info.to_excel(writer, index=False, sheet_name='all_info')
    writer.close()
    print(f"result was written to {output_path}")

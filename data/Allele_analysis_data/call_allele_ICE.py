# -*- coding: utf-8 -*-
# @Author  : xiang.zhang
# @FileName: call_allele_ICE.py
# @Software: PyCharm
# @usage   : call indels with ICE package
import warnings
warnings.filterwarnings("ignore")
import sys, os
project_path = r"D:\gitFolder\CasDelta"
sys.path.extend([project_path])
sys.path.extend([os.path.join(project_path, 'ice')])
from ice.ice.analysis import single_sanger_analysis
import os
import argparse
from tqdm import tqdm
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Analyze Sanger reads to infer crispr edit outcomes')
    parser.add_argument('--out', dest='out', help='Output directory path (defaults to .)', required=False, default=None)
    parser.add_argument('--data', dest='data', help='Data path, where .ab1 files are located (required)', required=True,
                        default=None)
    parser.add_argument('--wt', dest='wild type ab1', help = 'wild type ab1', required=True, default=None)
    parser.add_argument('--indel_max_size', dest='indel_max_size', help='max size of indel', default=50)
    parser.add_argument('--verbose', dest='verbose', action='store_true', help='Display verbose output', default=False)
    parser.add_argument('--allprops', dest='allprops', action='store_true', default=False,
                        help="Output all Edit Proposals, even if they have zero contribution")
    parser.add_argument('--target', dest='target',
                        help='Target sequence(s) (17-23 bases, RNA or DNA, comma separated), (REQUIRED)',
                        required=True)
    parser.add_argument('--donor', dest='donor', help='Donor DNA sequence for HDR (Optional)', required=False,
                        default=None)
    parser.add_argument('--target_name', dest='Name of target site', default='target1')
    args = parser.parse_args()

    ts_name = args.target_name
    out_path = r"E:\data_analysis\Projects\2024_08_29_Cas12_delta_Improve\write_paper\srcPaper\sangle_analysis_ICE"
    out_dir = os.path.join(out_path, ts_name)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    data_dir = args.data
    wt_ab1 = args.wt
    assert os.path.exists(data_dir) and os.path.exists(wt_ab1)

    indel_max_size = args.indel_max_size

    verbose = args.verbose
    allprops = args.allprops


    guide = args.target
    donor = args.donor

    files_as = [i for i in os.listdir(data_dir) if i.endswith('.ab1')]
    try:
        files_as.pop(files_as.index(os.path.basename(wt_ab1)))
    except ValueError:
        pass

    for idx, input_ab1 in enumerate(tqdm(files_as)):
        prefix = Path(input_ab1).stem
        single_sanger_analysis(control_path=wt_ab1,
                               sample_path=os.path.join(data_dir, input_ab1),
                               base_outputname=os.path.join(out_dir, prefix),
                               guide=guide,
                               donor=donor,
                               indel_max_size = indel_max_size,
                               verbose=verbose,
                               allprops=allprops)

if __name__ == '__main__':
    main()

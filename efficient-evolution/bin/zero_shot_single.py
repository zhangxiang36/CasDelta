# -*- coding: utf-8 -*-
# @Time    : 2025/10/17 17:21
# @Author  : xiang.zhang
# @FileName: zero_shot_single.py
# @Software: PyCharm
# @usage   : https://github.com/Gabriel-QIN/AI_Protein_Engneering/blob/master/zero_shot_single.py, calculating wt-marginal and masked marginal score

#!/usr/bin/python
import os
import os.path as osp
import argparse
import pathlib
import string
import torch
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import itertools
from typing import List, Tuple
from itertools import combinations
from transformers import EsmTokenizer, EsmForMaskedLM

# Credit: this code is adapted from https://github.com/amelie-iska/Variant-Effects.
# It contains several scripts for accessing variant effects (single-site and multi-site of mutations) with protein and DNA language models.
# To save computational cost, we only use wt-msarginals to infer variant effect in our study.

# Function to label a row in the DataFrame based on the scoring of mutations
def label_row(row, sequence, token_probs, tokenizer, offset_idx):
    # Extract wild type, index, and mutated type from the row
    wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
    assert sequence[idx] == wt, "The listed wildtype does not match the provided sequence"
    # Encode the wild type and mutated type
    wt_encoded, mt_encoded = tokenizer.encode(wt, add_special_tokens=False)[0], tokenizer.encode(mt, add_special_tokens=False)[0]
    # Calculate the score as the difference in log probabilities
    score = token_probs[0, 1 + idx, mt_encoded] - token_probs[0, 1 + idx, wt_encoded]
    return score.item()

# Function to compute pseudo-perplexity for a row, used in language model evaluation
def compute_pppl(row, sequence, model, tokenizer, offset_idx, device):
    wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
    assert sequence[idx] == wt, "The listed wildtype does not match the provided sequence"
    # Modify the sequence with the mutation
    sequence = sequence[:idx] + mt + sequence[(idx + 1) :]
    # Tokenize the modified sequence
    data = [("protein1", sequence)]
    batch_tokens = tokenizer.batch_encode_plus(data, return_tensors="pt", padding=True)["input_ids"]
    # Calculate log probabilities for each position
    log_probs = []
    for i in range(1, len(sequence) - 1):
        batch_tokens_masked = batch_tokens.clone()
        batch_tokens_masked[0, i] = tokenizer.mask_token_id
        with torch.no_grad():
            token_probs = torch.log_softmax(model(batch_tokens_masked.to(device)).logits, dim=-1)
        log_probs.append(token_probs[0, i, batch_tokens[0, i]].item())
    return sum(log_probs)

# Apply mutations on a given sequence
def apply_mutations(df, seq):
    seqlist = []
    for index,row in df.iterrows():
        pos = row['positions'].split('+')
        mut = row['mutant'].split('+')
        mutated_seq = list(seq)
        for p,aa in zip(pos, mut):
            mutated_seq[int(p)-1] = aa[-1]
        seqlist.append(''.join(mutated_seq))
    return seqlist

# Enumerate all of the single-site mutations given a protein sequence
def enumerate_mutation(seq, mutant_col='mutant'):
    data = []
    aalist = list("ACDEFGHIKLMNPQRSTVWY")
    # Enumerate whole sequence
    for seqidx in range(len(seq)):
        # Enumerate all of the 20 amino acids
        for aa in aalist:
            mt = aa[0]
            wt = seq[seqidx]
            # # Skip if the mutated residue is WT residue
            # if wt == mt:
            #     continue
            mutation = f'{wt}{seqidx+1}{mt}'
            data.append({
                f'{mutant_col}': mutation,
                'positions': f'{seqidx+1}',
            })
    df = pd.DataFrame(data)
    seqlist = apply_mutations(df, seq)
    # df['seq'] = seqlist
    return df

# Convert single-site mutant into a pandas DataFrame object
def get_mutated_seq(seq, mutation_list, mutant_col='mutant'):
    data = []
    for mutation in mutation_list:
        wt = mutation[0]
        mt = mutation[-1]

        position = int(mutation[1:-1])
        data.append({
            f'{mutant_col}': mutation,
            'positions': f'{position}',
        })
    df = pd.DataFrame(data)
    seqlist = apply_mutations(df, seq)
    # df['seq'] = seqlist
    return df

# Create and configure the argument parser for command line interface
def create_parser():
    parser = argparse.ArgumentParser(
        description="In silico deep mutational scan for single-site mutations using ESM predictions."
    )
    # Define the arguments the script can accept
    parser.add_argument('--name', '-n', type=str, default='example', help='Protein name.')
    parser.add_argument('--seq', '-s', type=str, required=True, help='Wild-type seq.')
    parser.add_argument('--mutation', '-mut', type=str, default=None, help='Mutation list')
    parser.add_argument('--plm', '-plm', type=str, default='esm2_t33_650M_UR50D', help='Name of the ESM model')
    parser.add_argument("--mutation-col", '-col', type=str, default="mutant", help="Column in the deep mutational scan labeling the mutation")
    parser.add_argument("--offset-idx", '-oi', type=int, default=1, help="Offset of the mutation positions")
    parser.add_argument("--scoring-strategy", '-ss', type=str, default="wt-marginals", choices=["wt-marginals", "masked-marginals", "pseudo-ppl"], help="Scoring strategy to use")
    parser.add_argument("--nogpu", action="store_true", help="Do not use GPU even if available")
    parser.add_argument('--gpu', '-gd', type=int, default=0, help='GPU device ID.')
    parser.add_argument('--topk', '-tk', type=int, default=20, help='Output top K mutations.')
    parser.add_argument('--output', '-o', type=str, default='output', help='Output directory for prediction results.')
    parser.add_argument('--prefix', '-pf', type=str, default='esm2', help='Output prefix.')
    return parser

# Main function to orchestrate mutation scoring process
def main(args):
    os.makedirs(f'{args.output}', exist_ok=True)
    # Determine to use GPU or CPU
    device = torch.device(f"cuda:{args.gpu}" if torch.cuda.is_available() and not args.nogpu else "cpu")
    # Load the chosen ESM model and tokenizer
    load_dir = args.plm
    model_name = args.plm
    model = EsmForMaskedLM.from_pretrained(model_name)
    for name, param in model.named_parameters():
        if 'contact_head.regression' in name:
            param.requires_grad = False
    tokenizer = EsmTokenizer.from_pretrained(model_name)
    model = model.to(device)
    model.eval()
    vocab = tokenizer.get_vocab()

    if args.mutation is not None:
        mutation_list = args.mutation.split(',')
        df = get_mutated_seq(args.seq, mutation_list, mutant_col=args.mutation_col)
    else:
        # enumerate all of the possible single-site mutations
        df = enumerate_mutation(args.seq, mutant_col=args.mutation_col)
    # Preprocess and encode the WT sequence
    data = [("protein1", args.seq)]
    batch_tokens = tokenizer.batch_encode_plus(data, return_tensors="pt", padding=True)["input_ids"].to(device)
    print(f'Generate {df.shape[0]} possible mutations.')
    # Apply selected scoring strategy
    if args.scoring_strategy == "wt-marginals":
        with torch.no_grad():
            token_probs = torch.log_softmax(model(batch_tokens).logits, dim=-1)
        df[args.plm] = df.apply(
            lambda row: label_row(row[args.mutation_col], args.seq, token_probs, tokenizer, args.offset_idx),
            axis=1
        )
    elif args.scoring_strategy == "masked-marginals":
        all_token_probs = []
        for i in tqdm(range(batch_tokens.size(1))):
            batch_tokens_masked = batch_tokens.clone()
            batch_tokens_masked[0, i] = tokenizer.mask_token_id
            with torch.no_grad():
                token_probs = torch.log_softmax(model(batch_tokens_masked).logits, dim=-1)
            all_token_probs.append(token_probs[:, i])
        token_probs = torch.cat(all_token_probs, dim=0).unsqueeze(0)
        df[args.plm] = df.apply(
            lambda row: label_row(row[args.mutation_col], args.seq, token_probs, tokenizer, args.offset_idx),
            axis=1
        )
    elif args.scoring_strategy == "pseudo-ppl":
        tqdm.pandas()
        df[args.plm] = df.progress_apply(
            lambda row: compute_pppl(row[args.mutation_col], args.seq, model, tokenizer, args.offset_idx, device),
            axis=1
        )
    df_top = df.sort_values(args.plm, ascending=False)
    topk = df_top[:args.topk]
    topk = topk.to_csv(osp.join(args.output, f'top{args.topk}_{args.name}.csv'), index=False, sep='\t')
    df[args.mutation_col] = [a[-1] for a in df[args.mutation_col].tolist()]
    df = df.pivot(index='positions', columns=args.mutation_col, values=args.plm)
    df['position'] = [int(a)+1 for a in range(len(df))]
    df['WT'] = [a for a in args.seq]
    df_sorted = df.sort_values(by='position')
    df_sorted = df_sorted[[df_sorted.columns[-2], df_sorted.columns[-1]] + list(df_sorted.columns[:-2])]
    # Save the scored mutations to a CSV file
    df_sorted.to_csv(osp.join(args.output, f'fitness_landscape_{args.prefix}_{args.name}.csv'), index=False, sep='\t')

# Script entry point for command-line interaction
if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)

# Example usage:
# myseq=MQRAVSVVARLGFRLQAFPPALCRPLSCAQEVLRRT && python zero_shot_single.py --name example  -s ${myseq} -o output/single

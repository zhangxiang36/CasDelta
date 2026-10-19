#!/usr/bin/env python3
"""
extract_CRISPRessoWGS_info.py

Scan CRISPRessoWGS output directory, extract key editing metrics from
quantification file and allele frequency table, and output a summary CSV.

Usage:
    python extract_CRISPRessoWGS_info.py \
        --input_dir /path/to/CRISPRessoWGS_output \
        --output summary_table.xlsx
"""

import argparse
import pandas as pd
from pathlib import Path

def parse_sample_region(dirname):
    """
    Try to extract sample/variant/target info from directory name.
    Customize this function based on your naming convention.

    Example dirnames:
        OsPi21_9_merge_crisprwgsout/CRISPRessoWGS_on_YG260025A005a_T436R

    Returns dict with keys: sample, variant, target, replicate (if applicable)
    """
    parts = dirname.parts
    target = '_'.join(parts[-2].split('_')[:2])
    event, variant = parts[-1].split('_')[-2:]
    # Default fallback: entire dirname as sample_id
    info = {
        'sample_id': event,
        'variant': variant,
        'target': target,
        'replicate': ''
    }

    return info

def read_quantification_file(filepath):
    """
    Read CRISPResso_quantification_of_editing_frequency.txt.

    Actual format (single-row, amplicon-level):
    Amplicon | Unmodified% | Modified% | Reads_in_input | Reads_aligned_all_amplicons |
    Reads_aligned | Unmodified | Modified | Discarded | Insertions | Deletions |
    Substitutions | Only Insertions | Only Deletions | Only Substitutions |
    Insertions and Deletions | Insertions and Substitutions |
    Deletions and Substitutions | Insertions Deletions and Substitutions
    """
    metrics = {}
    try:
        df = pd.read_csv(filepath, sep='\t')

        # Direct column mapping
        col_map = {
            'Unmodified%': 'unmodified_pct',
            'Modified%': 'modified_pct',
            'Reads_in_input': 'reads_input',
            'Reads_aligned': 'reads_aligned',
            'Reads_aligned_all_amplicons': 'reads_aligned_all',
            'Unmodified': 'reads_unmodified',
            'Modified': 'reads_modified',
            'Discarded': 'reads_discarded',
            'Insertions': 'reads_insertion',
            'Deletions': 'reads_deletion',
            'Substitutions': 'reads_substitution',
            'Only Insertions': 'reads_only_insertion',
            'Only Deletions': 'reads_only_deletion',
            'Only Substitutions': 'reads_only_substitution',
            'Insertions and Deletions': 'reads_ins_del',
            'Insertions and Substitutions': 'reads_ins_sub',
            'Deletions and Substitutions': 'reads_del_sub',
            'Insertions Deletions and Substitutions': 'reads_ins_del_sub'
        }

        for old_col, new_col in col_map.items():
            if old_col in df.columns:
                val = df[old_col].iloc[0]
                metrics[new_col] = val

        # Calculate derived metrics

        # 1. Indel-only frequency (excluding pure substitutions)
        only_del = metrics.get('reads_only_deletion', 0) or 0
        only_ins = metrics.get('reads_only_insertion', 0) or 0
        ins_del = metrics.get('reads_ins_del', 0) or 0
        total = metrics.get('reads_aligned', 0) or 1  # avoid div by zero

        indel_reads = only_del + only_ins + ins_del
        metrics['indel_read_count'] = indel_reads
        metrics['indel_read_pct'] = round(indel_reads / total * 100, 4)

        # 2. Pure substitution frequency
        only_sub = metrics.get('reads_only_substitution', 0) or 0
        metrics['substitution_only_pct'] = round(only_sub / total * 100, 4)

        # 3. Combined modification breakdown verification
        ins_sub = metrics.get('reads_ins_sub', 0) or 0
        del_sub = metrics.get('reads_del_sub', 0) or 0
        ins_del_sub = metrics.get('reads_ins_del_sub', 0) or 0

        total_modified = (only_ins + only_del + only_sub +
                          ins_del + ins_sub + del_sub + ins_del_sub)
        metrics['total_modified_reads_check'] = total_modified

        # 4. Insertion vs deletion ratio
        ins_total = only_ins + ins_del + ins_sub + ins_del_sub
        del_total = only_del + ins_del + del_sub + ins_del_sub
        metrics['insertion_total_reads'] = ins_total
        metrics['deletion_total_reads'] = del_total
        if ins_total + del_total > 0:
            metrics['ins_del_ratio'] = round(ins_total / (ins_total + del_total), 4)
        else:
            metrics['ins_del_ratio'] = None

        # 5. Modified% sanity check vs breakdown sum
        modified_pct = metrics.get('modified_pct', 0) or 0
        breakdown_pct = round(total_modified / total * 100, 4)
        metrics['modified_breakdown_pct'] = breakdown_pct
        metrics['modified_pct_diff'] = round(abs(modified_pct - breakdown_pct), 4)

    except Exception as e:
        print(f"  Warning: could not read {filepath}: {e}")

    return metrics

def read_allele_table(filepath, min_allele_pct=0.05):
    """
    Read Alleles_frequency_table.txt and extract allele-level metrics.

    Actual columns observed in your data:
    Aligned_Sequence, Reference_Sequence, Reference_Name, Read_Status,
    n_deleted, n_inserted, n_mutated, #Reads, %Reads

    Key observation: many MODIFIED rows have n_deleted=0, n_inserted=0,
    n_mutated>0 – these are PURE SUBSTITUTION events, not indels.
    """
    metrics = {}
    try:
        df = pd.read_csv(filepath, sep='\t')

        # Total reads represented in allele table
        if '#Reads' in df.columns:
            metrics['allele_table_reads'] = df['#Reads'].sum()

        # Filter out very low-frequency alleles (sequencing noise)
        df_filtered = df[df['%Reads'] >= min_allele_pct].copy()
        df = df_filtered

        # Separate unmodified and modified alleles
        unmod = df[df['Unedited']]
        mod = df[~df['Unedited']]

        # --- Basic percentages ---
        # Unmodified percentage
        if not unmod.empty:
            metrics['allele_unmodified_pct'] = unmod['%Reads'].sum()
        else:
            metrics['allele_unmodified_pct'] = 0.0

        # Modified (edited) percentage
        if not mod.empty:
            metrics['allele_modified_pct'] = mod['%Reads'].sum()
        else:
            metrics['allele_modified_pct'] = 0.0

        # Number of distinct edited alleles (after filtering)
        metrics['n_edited_alleles'] = len(mod)

        # --- Top edited allele ---
        if not mod.empty:
            metrics['top_edited_allele_pct'] = mod['%Reads'].max()
            # Get top allele details
            top_idx = mod['%Reads'].idxmax()
            top_allele = mod.loc[top_idx]
            metrics['top_allele_n_deleted'] = top_allele.get('n_deleted', None)
            metrics['top_allele_n_inserted'] = top_allele.get('n_inserted', None)
            metrics['top_allele_n_mutated'] = top_allele.get('n_mutated', None)
            # Truncate long sequences for readability
            seq = str(top_allele.get('Aligned_Sequence', ''))
            metrics['top_allele_seq'] = seq[:80] + '...' if len(seq) > 83 else seq
        else:
            metrics['top_edited_allele_pct'] = 0.0
            metrics['top_allele_n_deleted'] = None
            metrics['top_allele_n_inserted'] = None
            metrics['top_allele_n_mutated'] = None
            metrics['top_allele_seq'] = ''

        # --- Chimerism indicators ---
        if not mod.empty:
            top_frac = mod['%Reads'].max()

            # Dominant allele present (>80% of edited reads)
            metrics['has_dominant_allele'] = int(top_frac > 80)

            # Highly chimeric: many alleles, no dominant one
            metrics['is_highly_chimeric'] = int(
                (len(mod) >= 3) and (top_frac < 50)
            )

            # Biallelic-like: top 2 edited alleles account for >70% of edited reads
            top2 = mod.nlargest(2, '%Reads')['%Reads'].sum()
            metrics['biallelic_like'] = int((top2 > 70) and (len(mod) <= 3))
        else:
            metrics['has_dominant_allele'] = 0
            metrics['is_highly_chimeric'] = 0
            metrics['biallelic_like'] = 0

        # --- Indel vs substitution classification ---
        # CRITICAL: separate true indels from pure substitutions
        if not mod.empty:
            # True indels: n_deleted > 0 OR n_inserted > 0
            mod_with_indel = mod[(mod['n_deleted'] > 0) | (mod['n_inserted'] > 0)]
            # Pure substitutions: n_deleted=0, n_inserted=0, n_mutated>0
            mod_pure_sub = mod[(mod['n_deleted'] == 0) &
                               (mod['n_inserted'] == 0) &
                               (mod['n_mutated'] > 0)]
            # Complex/unclassified: other combinations (should be rare)
            mod_other = mod[~mod.index.isin(mod_with_indel.index) &
                            ~mod.index.isin(mod_pure_sub.index)]

            metrics['indel_allele_reads_pct'] = round(mod_with_indel['%Reads'].sum(),
                                                      4) if not mod_with_indel.empty else 0.0
            metrics['pure_substitution_allele_pct'] = round(mod_pure_sub['%Reads'].sum(),
                                                            4) if not mod_pure_sub.empty else 0.0
            metrics['other_mod_allele_pct'] = round(mod_other['%Reads'].sum(), 4) if not mod_other.empty else 0.0
            metrics['n_indel_alleles'] = len(mod_with_indel)
            metrics['n_pure_substitution_alleles'] = len(mod_pure_sub)
        else:
            metrics['indel_allele_reads_pct'] = 0.0
            metrics['pure_substitution_allele_pct'] = 0.0
            metrics['other_mod_allele_pct'] = 0.0
            metrics['n_indel_alleles'] = 0
            metrics['n_pure_substitution_alleles'] = 0

        # --- Frameshift calculation (only for true indels) ---
        if not mod_with_indel.empty:
            def is_frameshift(row):
                net = row.get('n_inserted', 0) - row.get('n_deleted', 0)
                return abs(net) % 3 != 0

            mod_with_indel = mod_with_indel.copy()
            mod_with_indel['frameshift'] = mod_with_indel.apply(is_frameshift, axis=1)
            fs_pct = mod_with_indel[mod_with_indel['frameshift']]['%Reads'].sum()
            metrics['frameshift_allele_pct'] = round(fs_pct, 4)
        else:
            metrics['frameshift_allele_pct'] = 0.0

        # --- Indel size statistics (from true indels only) ---
        if not mod_with_indel.empty:
            del_rows = mod_with_indel[mod_with_indel['n_deleted'] > 0]
            ins_rows = mod_with_indel[mod_with_indel['n_inserted'] > 0]

            metrics['mean_indel_size_deleted'] = round(del_rows['n_deleted'].mean(), 2) if not del_rows.empty else 0
            metrics['mean_indel_size_inserted'] = round(ins_rows['n_inserted'].mean(), 2) if not ins_rows.empty else 0
            metrics['max_indel_size_deleted'] = int(del_rows['n_deleted'].max()) if not del_rows.empty else 0
            metrics['max_indel_size_inserted'] = int(ins_rows['n_inserted'].max()) if not ins_rows.empty else 0
        else:
            metrics['mean_indel_size_deleted'] = 0
            metrics['mean_indel_size_inserted'] = 0
            metrics['max_indel_size_deleted'] = 0
            metrics['max_indel_size_inserted'] = 0

        # --- Small indel proportion (1-3 bp indels, often more reliable) ---
        if not mod_with_indel.empty:
            small_indel = mod_with_indel[
                ((mod_with_indel['n_deleted'] > 0) & (mod_with_indel['n_deleted'] <= 3)) |
                ((mod_with_indel['n_inserted'] > 0) & (mod_with_indel['n_inserted'] <= 3))
                ]
            metrics['small_indel_allele_pct'] = round(small_indel['%Reads'].sum(), 4) if not small_indel.empty else 0.0
        else:
            metrics['small_indel_allele_pct'] = 0.0
        if not unmod.empty:
            metrics['has_wt_allele'] = True
        else:
            metrics['has_wt_allele'] = False
    except Exception as e:
        print(f"  Warning: could not process {filepath}: {e}")

    return metrics

def get_qc_metrics(dir_path):
    """Get basic QC from mapping statistics.

    Actual format (tab-separated, single header + data row):
    READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\t...
    2294\t2294\t2294\t...
    """
    qc = {}
    map_file = dir_path / 'CRISPResso_mapping_statistics.txt'
    if map_file.exists():
        try:
            df = pd.read_csv(map_file, sep='\t')
            if not df.empty:
                row = df.iloc[0]
                # Map columns to our metric names
                col_map = {
                    'READS IN INPUTS': 'raw_reads',
                    'READS AFTER PREPROCESSING': 'clean_reads',
                    'READS ALIGNED': 'aligned_reads'
                }
                for orig_col, metric_name in col_map.items():
                    if orig_col in df.columns:
                        qc[metric_name] = int(row[orig_col])

                # Calculate alignment rate
                raw = qc.get('raw_reads', 0)
                aligned = qc.get('aligned_reads', 0)
                if raw > 0:
                    qc['alignment_rate'] = round(aligned / raw * 100, 2)
                else:
                    qc['alignment_rate'] = None

        except Exception as e:
            print(f"  Warning: could not read mapping stats: {e}")

    return qc


def main():
    parser = argparse.ArgumentParser(
        description='Extract CRISPRessoWGS summary table'
    )
    parser.add_argument(
        '--input_dir', '-i',
        required=True,
        help='Path to CRISPRessoWGS output directory'
    )
    parser.add_argument(
        '--output', '-o',
        default='crispressowgs_summary.xlsx',
        help='Output xlsx filename'
    )
    parser.add_argument(
        '--min_depth', type=int, default=200,
        help='Minimum aligned reads to include (default: 200)'
    )
    parser.add_argument(
        '--min_allele_pct', type=float, default=5,
        help='Minimum allele %% to keep in allele table (default: 5, means 5%)'
    )
    parser.add_argument(
        '--plant_lv_edited_pct', type=float, default=10.0,
        help='Plant-level edited allele percentage threshold (default: 10)'
    )

    args = parser.parse_args()

    input_path = Path(args.input_dir)
    if not input_path.is_dir():
        print(f"Error: {args.input_dir} is not a valid directory")
        return

    results = []
    skipped = []

    print(f"Scanning {input_path} ...")

    # Iterate through subdirectories (each is a CRISPRessoWGS region/sample)
    for subdir in sorted(input_path.iterdir()):
        if not subdir.is_dir():
            continue

        dirname = subdir.name
        info_dir = [i for i in subdir.iterdir() if i.is_dir() and i.parts[-1].startswith('CRISPResso')][0]
        print(f"\nProcessing: {dirname}")

        # Parse sample/variant/target from directory name
        info = parse_sample_region(subdir)

        # Check for required files
        quant_file = info_dir / 'CRISPResso_quantification_of_editing_frequency.txt'
        assert sorted(info_dir.glob('Alleles_frequency_table_around_sgRNA_*.txt')).__len__() == 1
        allele_file = info_dir / sorted(info_dir.glob('Alleles_frequency_table_around_sgRNA_*.txt'))[0]

        if not quant_file.exists():
            print(f"  Skipping: no quantification file found")
            skipped.append({'directory': dirname, 'reason': 'no quantification file'})
            continue

        if not allele_file.exists():
            print(f"  Skipping: no allele frequency table found")
            skipped.append({'directory': dirname, 'reason': 'no allele table'})
            continue
        print('find info files')
        # Extract metrics
        quant_metrics = read_quantification_file(quant_file)
        allele_metrics = read_allele_table(allele_file, args.min_allele_pct)
        qc_metrics = get_qc_metrics(info_dir)

        # Combine all metrics
        row = {
            'directory': dirname,
            'sample_id': info['sample_id'],
            'variant': info['variant'],
            'target': info['target'],
            'replicate': info.get('replicate', ''),
        }
        row.update(qc_metrics)
        row.update(quant_metrics)
        row.update(allele_metrics)

        # Add plant-level edit call (for compatibility with existing data)
        if 'indel_read_pct' in row and row['indel_read_pct'] is not None:
            row['plant_level_edited'] = int(row['indel_read_pct'] > args.plant_lv_edited_pct)
        else:
            row['plant_level_edited'] = None

        # Check depth filter
        aligned = row.get('reads_aligned', 0)
        if aligned and aligned < args.min_depth:
            print(f"  Warning: low depth ({aligned} aligned reads), including but flagging")
            row['low_depth_flag'] = 1
        else:
            row['low_depth_flag'] = 0

        results.append(row)
        print(f"  OK: indel_read_pct={row.get('indel_read_pct', 'NA'):.2f}%, "
              f"modified_pct={row.get('modified_pct', 'NA'):.2f}%, "
              f"indel_allele_pct={row.get('indel_allele_reads_pct', 'NA'):.2f}%, "
              f"sub_allele_pct={row.get('pure_substitution_allele_pct', 'NA'):.2f}%, "
              f"top_allele={row.get('top_edited_allele_pct', 'NA'):.1f}%, "
              f"n_alleles={row.get('n_edited_alleles', 'NA')}")

    # Create DataFrame
    if not results:
        print("\nNo data extracted. Check directory structure.")
        return

    df = pd.DataFrame(results)

    # Define column order
    base_cols = ['directory', 'sample_id', 'variant', 'target', 'replicate']
    qc_cols = ['raw_reads', 'clean_reads', 'aligned_reads', 'alignment_rate']
    quant_cols = [
        'reads_input', 'reads_aligned', 'reads_aligned_all',
        'reads_unmodified', 'reads_modified', 'reads_discarded',
        'unmodified_pct', 'modified_pct',
        'indel_read_count', 'indel_read_pct',
        'substitution_only_pct', 'modified_breakdown_pct', 'modified_pct_diff',
        'reads_insertion', 'reads_deletion', 'reads_substitution',
        'reads_only_insertion', 'reads_only_deletion', 'reads_only_substitution',
        'reads_ins_del', 'reads_ins_sub', 'reads_del_sub', 'reads_ins_del_sub',
        'insertion_total_reads', 'deletion_total_reads', 'ins_del_ratio',
    ]
    allele_cols = [
        'allele_table_reads',
        'allele_unmodified_pct', 'allele_modified_pct',
        'n_edited_alleles', 'top_edited_allele_pct',
        'top_allele_n_deleted', 'top_allele_n_inserted', 'top_allele_n_mutated',
        'top_allele_seq',
        'has_dominant_allele', 'has_wt_allele', 'is_highly_chimeric', 'biallelic_like',
        'indel_allele_reads_pct', 'pure_substitution_allele_pct', 'other_mod_allele_pct',
        'n_indel_alleles', 'n_pure_substitution_alleles',
        'frameshift_allele_pct',
        'mean_indel_size_deleted', 'mean_indel_size_inserted',
        'max_indel_size_deleted', 'max_indel_size_inserted',
        'small_indel_allele_pct'
    ]
    extra_cols = ['plant_level_edited', 'low_depth_flag']

    # Keep only columns that exist in the dataframe
    ordered_cols = [c for c in base_cols + qc_cols + quant_cols + allele_cols + extra_cols
                    if c in df.columns]

    # Add any additional columns not in the predefined list
    extra = [c for c in df.columns if c not in ordered_cols]
    ordered_cols += extra

    df = df[ordered_cols]

    # Sort by variant, target, sample_id
    sort_cols = [c for c in ['variant', 'target', 'sample_id'] if c in df.columns]
    if sort_cols:
        df = df.sort_values(sort_cols)

    # Write output
    select_col = [
        "sample_id", "variant", "target", "replicate",
        "aligned_reads", "unmodified_pct", "modified_pct", "indel_read_count", "indel_read_pct", "n_edited_alleles",
        "indel_allele_reads_pct", "n_indel_alleles", "has_wt_allele","top_edited_allele_pct", "has_dominant_allele",
        "is_highly_chimeric", "biallelic_like", "mean_indel_size_deleted",
        "mean_indel_size_inserted", "max_indel_size_deleted", "max_indel_size_inserted", "plant_level_edited",
        "low_depth_flag"
    ]
    writer = pd.ExcelWriter(args.output)
    df[select_col].to_excel(writer, index=False, sheet_name='information')
    writer.close()
    print(f"\nDone! Summary written to: {args.output}")
    print(f"Total regions processed: {len(df)}")
    print(f"Skipped directories: {len(skipped)}")

    if skipped:
        skip_df = pd.DataFrame(skipped)
        skip_file = args.output.replace('.xlsx', '_skipped.xlsx')
        skip_df.to_excel(skip_file, index=False)
        print(f"Skipped directories logged to: {skip_file}")


if __name__ == '__main__':
    main()

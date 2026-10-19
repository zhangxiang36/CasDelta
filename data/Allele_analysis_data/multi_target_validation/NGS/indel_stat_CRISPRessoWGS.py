# -*- coding: utf-8 -*-
"""
统计 CRISPResso 等位基因频率表中各 deletion size 的占比。
规则：
  1. 用 n_deleted 作为 del_len（相比野生型的 deletion size）
  2. Unedited=True 的行视为野生型
  3. 去除占比 (%Reads) 低于阈值的条目（默认 5%，可调）
  4. 同 del_len 同 Unedited 状态求和
输出：
  a) del_len-是否为野生型-占比   (分号分隔)
  b) del_len,1;-51;-18;-17 形式 (野生型用1，突变用负del值)
"""
import argparse
import sys
from collections import defaultdict


def classify(n_types, has_wt):
    """按类型数和是否含 WT 判定倍型。"""
    if n_types >= 3:
        return "chimeric"
    if n_types == 2:
        return "het" if has_wt else "bi-allelic"
    # n_types == 1
    if n_types == 1 and not has_wt:
        return 'homo'
    else:
        raise ValueError('error in classify')

def load_data(path):
    """
    返回 [(net_len, unedited, pct), ...]
    真实净改变长度 = n_deleted - n_inserted
      > 0  → 净缺失，记为负数（缺失突变）
      < 0  → 净插入，记为正数（插入突变）
      = 0  → 既无缺失也无插入（纯 SNP / 完全匹配 / 野生型）
    """
    rows = []
    with open(path, encoding="utf-8") as f:
        header = f.readline()
        cols = header.rstrip("\n").split("\t")
        # 列顺序：Aligned_Sequence Reference_Sequence Unedited n_deleted n_inserted n_mutated #Reads %Reads
        idx = {c: i for i, c in enumerate(cols)}
        for line in f:
            if not line.strip():
                continue
            fld = line.rstrip("\n").split("\t")
            try:
                unedited = fld[idx["Unedited"]].strip().lower() == "true"
                n_del = int(fld[idx["n_deleted"]].strip())
                n_ins = int(fld[idx["n_inserted"]].strip())
                pct = float(fld[idx["%Reads"]].strip())
            except (ValueError, KeyError, IndexError) as e:
                print(f"[警告] 跳过无法解析的行: {line.rstrip()}: {e}", file=sys.stderr)
                continue
            # 净改变长度：缺失为正、插入为负
            net_change = n_del - n_ins
            # 最终编码：插入突变→正数，缺失突变→负数，无indel→0
            net_len = -net_change
            rows.append((net_len, unedited, pct, n_del, n_ins))
    return rows


def judgeAllele(crispresso_tsv, threshold = 5, wt_as_class = True):
    """

    :param crispresso_tsv: Alleles_frequency_table_around_sgRNA_XXXXX.txt
    :param threshold: exclude allele if it's percentage is lower than threshold
    :param wt_as_class: WT allele as a class
    :return: fmt1 (del_len-is WT?-percentage_of_allele), fmt2 (del_len1;del_len2), has_wt, genotype (het, homo, bi-allelic, chemiric)
    """
    rows = load_data(crispresso_tsv)

    # 按 net_len + Unedited 归并，过滤低于阈值的条目
    bucket = defaultdict(float)
    for net_len, unedited, pct, n_del, n_ins in rows:
        if pct < threshold:
            continue
        bucket[(net_len, unedited)] += pct

    if not bucket:
        print(f"[警告] 过滤阈值 {threshold}% 下没有数据，请调低阈值。", file=sys.stderr)
        sys.exit(0)

    # 排序：野生型排最前，再按净长度升序（负值缺失在前，正值插入在后）
    def sort_key(k):
        nl, ue = k
        return (0 if ue else 1, nl)
    items = sorted(bucket.items(), key=lambda kv: sort_key(kv[0]))

    # 形式一：net_len-是否为野生型-占比
    #   net_len 为正数=插入突变，负数=缺失突变，0=无indel/野生型
    fmt1 = ";".join(f"{nl}-{str(ue)}-{p:.2f}" for (nl, ue), p in items)

    # 形式二：野生型用 1，插入突变用正数，缺失突变用负数
    parts = []
    for (nl, ue), p in items:
        if ue:
            continue
        else:
            parts.append(f"{nl}")
    fmt2 = ";".join(parts)

    print("=== 形式一：net_len-是否为野生型-占比 ===")
    print("（net_len = n_deleted - n_inserted；正数=插入突变，负数=缺失突变，0=无indel/野生型）")
    print(fmt1)
    print()
    print("=== 形式二：del_len,1;-51;-18;-17 ===")
    print("（野生型→1，插入突变→正数，缺失突变→负数）")
    print(fmt2)
    print()

    # 明细表格
    print("=== 明细（归并后）===")
    print(f"{'net_len':>8} {'类型':>10} {'野生型':>6} {'占比(%)':>10}")
    for (nl, ue), p in items:
        if ue:
            label = "WT/无indel"
        elif nl > 0:
            label = f"插入+{nl}"
        elif nl < 0:
            label = f"缺失{nl}"
        else:
            label = "无indel突变"
        print(f"{nl:>8} {label:>10} {int(ue):>6} {p:>10.4f}")

        # ===== 倍型分类 =====
        # 关键：先按"是不是 WT"分两大类，再在突变里按 net_len 细分
    mut_netlens = set()
    wt_total = 0.0
    mut_total = 0.0
    for (nl, ue), p in items:
        if ue:
            wt_total += p
        else:
            mut_netlens.add(nl)
            mut_total += p

    has_wt = wt_total > 0
    n_mut_alleles = len(mut_netlens)

    # 类型数：默认 WT 算一种；可关闭
    if wt_as_class:
        n_types = (1 if has_wt else 0) + n_mut_alleles
    else:
        n_types = n_mut_alleles

    genotype = classify(n_types, has_wt)
    return fmt1, fmt2, has_wt, genotype

if __name__=="__main__":

    p = argparse.ArgumentParser(description="统计 deletion size 占比")
    p.add_argument("tsv", help="CRISPResso 等位基因频率表（tab 分隔）")
    p.add_argument("-t", "--threshold", type=float, default=5.0,
                   help="最低 %%Reads 阈值，低于此值的条目被过滤（默认 5）")
    p.add_argument("--wt_as_class", action="store_true",
                   help="野生型作为一个allele类别进行统计")

    args = p.parse_args()

    crispresso_tsv = args.tsv
    threshold = args.threshold
    wt_as_class = args.wt_as_class

    result = judgeAllele(crispresso_tsv, threshold, wt_as_class)
    print(crispresso_tsv, result)

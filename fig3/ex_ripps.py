# -*- coding: utf-8 -*-
# @Time :2025/2/8 13:37
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Extract ripps sequence from ".gbk" files predicted by antiSMASH 7.0
The gbk files exceeds the maximum file limit allowed by GitHub
so they can be find at https://zenodo.org/uploads/16420027
unzip the compressed file and the following five non-metagenome gbk files are used to extract ripp sequences
GCA_044360605.1_ASM4436060v1_genomic.gbk
GCA_044360685.1_ASM4436068v1_genomic.gbk
GCA_000335475.2_ASM33547v2_genomic.gbk
GCA_002305875.1_ASM230587v1_genomic.gbk
GCA_001887355.1_ASM188735v1_genomic.gbk
'''
import sys
sys.path.append('../')
import os
import pandas as pd
from astool.antismash_utils import AntismashRegionGBKParser

if __name__ == '__main__':
    tsv_output = 'ripps_cystobacter.tsv'
    gbk_files = os.listdir('amp_gbk')
    gbk_input = ['./amp_gbk/' + file for file in gbk_files]


    gbk_name_ls = []
    leader_seqs = []
    core_seqs = []
    gbk_files = []
    for gbk_dir in gbk_input:
        gbk_file = AntismashRegionGBKParser(gbk_dir)
        leader_seq, core_seq = gbk_file.ex_lanthipeptide()
        for ls, cs in zip(leader_seq, core_seq):
            gbk_name_ls.append(os.path.basename(gbk_dir))
            leader_seqs.append(ls)
            core_seqs.append(cs)
            gbk_files.append(gbk_dir)

    pd.DataFrame({
        "gbk_file": gbk_files,
        "gbk_name": gbk_name_ls,
        "leader_seq": leader_seqs,
        "core_seq": core_seqs
    }).to_csv(tsv_output, sep='\t', index=False)


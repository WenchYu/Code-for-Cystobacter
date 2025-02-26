# -*- coding: utf-8 -*-
# @Time :2025/2/8 13:37
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import os
from script import extractGBK
import pandas as pd
from astool.antismash_utils import AntismashRegionGBKParser
from astool.utils import get_gbk_dir_ls

if __name__ == '__main__':
    tsv_output = 'ripps_cystobacter.tsv'
    gbk_files = os.listdir('gbk')
    gbk_input = ['./gbk/' + file for file in gbk_files]

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


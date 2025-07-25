# -*- coding: utf-8 -*-
# @Time :2025/5/5 16:53
# @Auther :Yuwenchao
# @Software : PyCharm
'''
header of each sequence must be pure number (not string)
Conver to fasta format
'''
import time

import pandas as pd

if __name__ == '__main__':
    t = time.time()

    rc_csv = 'ripps_cystobacter.tsv'
    rc_fasta = 'ripps_cystobacter.fasta'

    rc_df = pd.read_csv(rc_csv,sep='\t') # 'gbk_name', 'amp_prediction_index'
    print(rc_df.gbk_name)
    print(rc_df.amp_prediction_index)

    with open(rc_fasta,'r+') as f:
        text = f.readline()

    print(f'Finished in {(time.time() - t) / 60:.2f} min')

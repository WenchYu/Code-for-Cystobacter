# -*- coding: utf-8 -*-
# @Time :2025/2/8 16:09
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''

import pandas as pd

if __name__ == '__main__':
    aa_list =['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    replace = 'X'
    tsv_file = 'ripps_cystobacter.tsv'
    df = pd.read_csv(tsv_file,sep='\t')
    gbk_names = df.gbk_name.tolist()
    core_seqs = df.core_seq.tolist()
    with open('ripps_cystobacter.fasta', 'w') as f:
        for i in range(len(gbk_names)):
            core_seqs[i] = core_seqs[i].replace(' ', '')
            core_seqs[i] = ''.join([char if char in aa_list else 'X' for char in core_seqs[i]])
            f.write(f'>{gbk_names[i]}\n{core_seqs[i]}\n')





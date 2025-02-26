# -*- coding: utf-8 -*-
# @Time :2025/2/12 16:41
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import collections

import pandas as pd

from my_packages import functions

if __name__ == '__main__':
    file = 'Fig2A_BGC-SSN_attributes.csv'
    df = functions.df_preprocess(file)
    clusters = list(collections.Counter(df.__ccCluster).keys())
    bigscape_class = list(collections.Counter(df['BiG-SCAPE class']).keys())

    counts_df = pd.DataFrame(columns=['NRPS', 'Terpene', 'RiPPs', 'Others', 'PKS-NRP_Hybrids', 'Saccharides', 'PKSother', 'PKSI'])
    counts_df.loc[0] = [0] * 8

    for cluster in clusters:
        cluster_df = df[df.__ccCluster == cluster]
        num_of_strains = len(collections.Counter(cluster_df.Strains))
        cla =  cluster_df['BiG-SCAPE class'].tolist()[0]
        if num_of_strains == 1:
            counts_df.loc[0, cla] +=1
    print(counts_df)


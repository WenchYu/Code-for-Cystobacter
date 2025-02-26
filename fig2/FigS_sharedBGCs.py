# -*- coding: utf-8 -*-
# @Time :2025/2/12 15:49
# @Auther :Yuwenchao
# @Software : PyCharm
'''
hist
'''
import collections
from my_packages import functions


if __name__ == '__main__':
    file = 'Fig2A_BGC-SSN_attributes.csv'
    df = functions.df_preprocess(file)
    # strains = list(collections.Counter(df.Strains).keys())
    clusters = list(collections.Counter(df.__ccCluster).keys())
    bigscape_class = list(collections.Counter(df['BiG-SCAPE class']).keys())
    five, four, three, two, one = 0,0,0,0,0
    print(len(clusters))
    for cluster in clusters:
        cluster_df = df[df.__ccCluster == cluster]
        num_of_strains = len(collections.Counter(cluster_df.Strains))
        if num_of_strains >= 5: five +=1
        if num_of_strains == 4: four +=1
        if num_of_strains == 3: three +=1
        if num_of_strains == 2: two +=1
        if num_of_strains == 1: one +=1

    print(five, four, three, two, one)


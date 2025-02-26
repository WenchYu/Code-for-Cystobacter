# -*- coding: utf-8 -*-
# @Time :2025/2/9 21:34
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import time
from my_packages import functions

if __name__ == '__main__':

    file = 'Fig2B_MCS-SSN_activity.xlsx'
    df = functions.df_preprocess(file)

    '''New columns'''
    new_columns = []
    for i in range(len(df)):
        try:
            new_columns.extend(df.loc[i,'activity'].split('/'))
        except:
            new_columns.append(df.loc[i,'activity'])
    new_columns = list(set(new_columns))
    for col in new_columns:
        df[col] = 0

    for i in range(len(df)):
        try:
            list = (df.loc[i,'activity'].split('/'))
            for element in list:
                df.loc[i,element] = df.loc[i,element]+1
        except:
            df.loc[i, df.loc[i,'activity']] = df.loc[i, df.loc[i,'activity']] + 1

    df.to_csv('Fig2B_MCS-SN_NodeAttributes.csv')



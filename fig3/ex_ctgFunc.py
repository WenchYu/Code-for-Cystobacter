# -*- coding: utf-8 -*-
# @Time :2025/2/11 14:10
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import os
import json
import pprint
import pandas as pd


def ex_ctgFunc(file, ):
    '''
    有些的key不对
    :param file:
    :return:
    '''
    with open(file, 'r') as f:
        data = json.load(f)
        # data.keys: ['version', 'input_file', 'records', 'timings', 'taxon', 'schema']

    level1 = data['records'][0]
    # dict_keys ['id', 'seq', 'features', 'name', 'description', 'dbxrefs', 'annotations', 'letter_annotations', 'areas', 'gc_content', 'modules']

    level2 = level1['modules']  # list containing 12956 elements
    # dict_keys =
    # (['antismash.detection.hmm_detection', 'antismash.detection.cluster_hmmer', 'antismash.detection.genefunctions', 'antismash.detection.nrps_pks_domains', 'antismash.detection.tigrfam', 'antismash.modules.active_site_finder', 'antismash.modules.cluster_compare', 'antismash.modules.clusterblast', 'antismash.modules.lanthipeptides', 'antismash.modules.lassopeptides', 'antismash.modules.nrps_pks', 'antismash.modules.pfam2go', 'antismash.modules.rrefinder', 'antismash.modules.sactipeptides', 'antismash.modules.t2pks', 'antismash.modules.tfbs_finder', 'antismash.modules.thiopeptides', 'antismash.modules.tta'])

    level3 = level2['antismash.detection.cluster_hmmer']
    # dict_keys(['hits', 'record id', 'schema', 'max evalue', 'min score', 'database', 'tool'])

    locuss, descriptions = [], []
    level4 = level3['hits']  # list[{},{},{}]
    for i in range(len(level4)):
        dict = level4[i]
        locus = dict['locus_tag']
        locuss.append(locus)
        description = dict['description']
        descriptions.append(description)
    df = pd.DataFrame({'locus_tag': locuss, 'description': descriptions})
    df.to_csv(f'ctgFunc{os.path.basename(file)}.csv')
    print('Done!!!')

if __name__ == '__main__':

    files = [os.path.join('./total_gbk/',file) for file in os.listdir('./total_gbk')]
    for file in files:
        try:
            ex_ctgFunc(file)
        except:print(file)

    print(f'')

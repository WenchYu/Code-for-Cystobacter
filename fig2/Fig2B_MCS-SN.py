# -*- coding: utf-8 -*-
# @Time :2025/2/9 00:22
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import sys
sys.path.append('../my_packages')
import pandas as pd
import networkx as nx
from rdkit.Chem import rdFMCS,MolFromSmiles
from tqdm import trange
from my_packages import functions

def df_preprocess(filename):
    '''
    Preprocess DataFrame by removing empty columns and resetting index.
    '''
    if filename.endswith('.csv'):
        df = pd.read_csv(filename, low_memory=False)
    elif filename.endswith('.tsv'):
        df = pd.read_csv(filename, sep='\t', low_memory=False)
    elif filename.endswith('.xlsx') or filename.endswith('.xls'):
        df = pd.read_excel(filename)
    else:
        raise ValueError("Unsupported file format. Please use .csv, .tsv, or .xlsx files.")

    if  df.index[-1] != len(df)-1:
        df.index.name = ''
        df.reset_index(inplace=True)
    return df

def MCS(mol1, mol2):

    mcs = rdFMCS.FindMCS([mol1, mol2]
                         , bondCompare=rdFMCS.BondCompare.CompareOrder
                         , atomCompare=rdFMCS.AtomCompare.CompareAny
                         , maximizeBonds = False
                         , ringMatchesRingOnly=False
                         , matchValences=False
                         , timeout=10
                         )

    mcs_num_bonds = mcs.numBonds
    mol1_num_bonds = mol1.GetNumBonds()
    mol2_num_bonds = mol2.GetNumBonds()

    similarity = mcs_num_bonds / ((mol1_num_bonds + mol2_num_bonds) - mcs_num_bonds)
    return similarity

if __name__ == '__main__':
    file = './/cystobacter.xlsx'

    threshold = 0.7
    topk = 10
    output = f"MCS-SN_Sim{threshold}Topk{topk}.graphml"
    df = df_preprocess(file)

    ids = df.npaid.tolist()

    names = df.compound_name.tolist()
    smiles = df.compound_smiles.tolist()
    bgcs = df.mibig_ids.tolist()
    strains = df.strain
    compound_level1s = df['Parent Level 1'].tolist()

    G = nx.Graph()
    for i in range(len(ids)):  # create nodes
        node_attr = {'names': names[i],
                     'smile': smiles[i],
                     'bgc':bgcs[i],
                     'strain':strains[i],
                     'parent_level1': compound_level1s[i]
                     }
        G.add_node(ids[i], **node_attr)

    for i in trange(len(ids)):
        smile1 = smiles[i]
        mol1 = MolFromSmiles(smile1)

        for j in range(i-1):
            smile2 = smiles[j]
            mol2 = MolFromSmiles(smile2)

            sim = MCS(mol1,mol2)
            if sim >= threshold:
                edge_attr = {'pair_similarity': sim}
                G.add_edge(ids[i], ids[j], **edge_attr)
        G = functions.mn_curating(G, topk)
        nx.write_graphml(G, output)

    print('')
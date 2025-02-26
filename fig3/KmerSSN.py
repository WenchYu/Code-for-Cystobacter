

import networkx as nx
from my_packages import functions

def parse_polymer_file(inp) -> dict:
    """Parse the input polymer file and convert data to dictionary format."""
    with open(inp, 'r') as f:
        data_ls = f.readlines()
    polymer_dict = {}
    for idx in range(len(data_ls)):
        if data_ls[idx].startswith('>'):
            header = data_ls[idx].strip()[1:]
            polymer_ls = data_ls[idx + 1].strip().split('-')
            polymer_dict[header] = polymer_ls
    return polymer_dict

def calculate_kmer_frequency(sequence, k):
    """
    计算给定序列中的k-mer频率。

    参数:
        sequence (str): 输入的序列字符串。
        k (int): k-mer的长度。

    返回:
        dict: 包含k-mer及其频率的字典。
    """
    if k <= 0:
        raise ValueError("k必须大于0")
    n = len(sequence)
    if n < k:
        return {}  # 序列长度小于k，无法生成k-mer
    kmers = [sequence[i:i + k] for i in range(n - k + 1)]
    kmer_counts = {}
    for kmer in kmers:
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1
    total_kmers = len(kmers)
    kmer_frequency = {kmer: count / total_kmers for kmer, count in kmer_counts.items()}
    return kmer_frequency

def calculate_weighted_euclidean_distance(freq1, freq2, weights=None):
    """
    计算两个k-mer频率字典之间的加权欧氏距离。

    参数:
        freq1 (dict): 第一个k-mer频率字典。
        freq2 (dict): 第二个k-mer频率字典。
        weights (dict, optional): 每个k-mer的权重字典。默认为None。

    返回:
        float: 两个频率字典之间的加权欧氏距离。
    """
    # 获取所有唯一的k-mer
    all_kmers = set(freq1.keys()).union(freq2.keys())

    # 如果未提供权重，则默认为1
    if weights is None:
        weights = {kmer: 1 for kmer in all_kmers}

    # 计算加权欧氏距离
    distance = 0.0
    for kmer in all_kmers:
        f1 = freq1.get(kmer, 0.0)
        f2 = freq2.get(kmer, 0.0)
        weight = weights.get(kmer, 1.0)  # 如果kmer不在权重中，默认权重为1
        distance += ((f1 - f2) * weight) ** 2
    distance = distance ** 0.5  # 欧氏距离需要开平方
    return 1-distance



def clustering(polymer_dict, k, threshold, output, topk):
    keys = []
    values = []
    for key, value in polymer_dict.items():
        keys.append(key)
        values.extend(value)

    G = nx.Graph()
    for i in range(len(keys)):  # create nodes
        node_attr = {'gbk_name': keys[i], 'core_peptide': values[i]}
        G.add_node(str(i), **node_attr)

    for i in range(len(keys)):  # create edges
        seq1 = values[i]
        freq1 = calculate_kmer_frequency(seq1, k)

        for j in range(i - 1):
            seq2 = values[j]
            freq2 = calculate_kmer_frequency(seq2, k)
            sim = calculate_weighted_euclidean_distance(freq1, freq2)

            if sim >= threshold:
                edge_attr = {'pair_similarity': sim}
                G.add_edge(str(i), str(j), **edge_attr)
    G = functions.mn_curating(G, topk)
    nx.write_graphml(G, output)

if __name__ == '__main__':

    w = 1
    threshold = 0.9
    k = 2
    topk = 5
    output = f'CorePeptide_K{k}Sim{threshold}Top{topk}.graphml'
    inp = 'ripps_cystobacter.fasta'
    polymer_dict = parse_polymer_file(inp)
    clustering(polymer_dict,k,threshold,output,topk)





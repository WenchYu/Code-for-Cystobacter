
import sys
sys.path.append('../')
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
    Inject language or reference

    param:
    sequence (str): input string
    k (int): length of k-mer

    returen:
        dict: dictionary of k-mer and frequency
    """
    if k <= 0:
        raise ValueError("k must > 0")
    n = len(sequence)
    if n < k:
        return {}
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
    Compute weighted Euclidean distance between two k-mer frequency dictionaries

    :param:
        freq1 (dict): The first k-mer frequency dictionary
        freq2 (dict): The second k-mer frequency dictionary.
        weights (dict, optional): the weight of each frequency, default = 1。

    :return:
        float: similarity (1-distance)
    """

    # Get all unique k-mers
    all_kmers = set(freq1.keys()).union(freq2.keys())


    if weights is None:
        weights = {kmer: 1 for kmer in all_kmers} # If no weight is provided, it defaults to 1

    # calculate distance
    distance = 0.0
    for kmer in all_kmers:
        f1 = freq1.get(kmer, 0.0)
        f2 = freq2.get(kmer, 0.0)
        weight = weights.get(kmer, 1.0)
        distance += ((f1 - f2) * weight) ** 2
    distance = distance ** 0.5
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
    threshold = 0.8
    k = 2
    topk = 5
    output = f'CorePeptide_K{k}Sim{threshold}Top{topk}.graphml'
    inp = 'ripps_cystobacter.fasta'
    polymer_dict = parse_polymer_file(inp)
    clustering(polymer_dict,k,threshold,output,topk)





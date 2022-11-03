# Plotting functions for project
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
from analysis_utils import call_aa_simple


def plot_obs_fitness_heatmap(reference_seq, double_mut_seq_list, obs_fitness_list):
    """
    Plots a heatmap of observed fitness on a len(reference_seq) x len(reference_seq) grid

    :param reference_seq: string of amino acids of reference protein sequence
    :param double_mut_seq_list: list of sequences of mutants with double mutations
    :param obs_fitness_list: list of observed fitness scores
    :return: None
    """
    obs_fitness_heatmap = np.ones((len(reference_seq), len(reference_seq))) * -6

    obs_fitness_pos_score_list = []

    for double_mut in range(len(double_mut_seq_list)):
        obs_fitness_pos_score_list_element = []
        sequence = double_mut_seq_list[double_mut]
        _, pos, _ = call_aa_simple(reference_seq, sequence)
        obs_fitness = obs_fitness_list[double_mut]

        # Append list element with double mutation positions and corresponding epistatic score / observed fitness
        obs_fitness_pos_score_list_element.append(pos)
        obs_fitness_pos_score_list_element.append(obs_fitness)

        # Append total epistatic score list with element
        obs_fitness_pos_score_list.append(obs_fitness_pos_score_list_element)

    for double_mut in range(len(obs_fitness_pos_score_list)):
        index_pair = obs_fitness_pos_score_list[double_mut][0]
        score_of_index_pair = obs_fitness_pos_score_list[double_mut][1]
        obs_fitness_heatmap[index_pair[0] - 1, index_pair[1] - 1] = score_of_index_pair
        obs_fitness_heatmap[index_pair[1] - 1, index_pair[0] - 1] = score_of_index_pair

    plt.imshow(obs_fitness_heatmap, cmap='viridis', interpolation='nearest')
    plt.colorbar()
    plt.show()


def plot_node_degree_distribution(epistasis_graph):
    """
    Given an epistasis graph, plot the node degree distribution (# epistatic interaction per amino acid position) for
    each node (amino acid position)

    :param epistasis_graph: Networkx graph
    :return: None
    """

    node_degree_list = list(map(list, sorted(epistasis_graph.degree, key=lambda x: x[1], reverse=True)))

    # Plot node degree list as bar chart
    plt.bar(np.array(node_degree_list)[:, 0], np.array(node_degree_list)[:, 1])
    plt.title("Node Degree Distribution")
    plt.xlabel("Amino Acid Position")
    plt.ylabel("Number of Epistatic Interactions")
    plt.show()

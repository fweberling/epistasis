import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import networkx as nx
from analysis_utils import preprocessing, double_mut_pos, epistatic_interaction_double_mutation, epistasis_graph, \
    call_aa_simple

# Specify order of mutations to be analysed
num_mut = 2

# Upload input files into panda data frame
data_frame = pd.read_csv('CPA_merge_filtered.csv')

# Specify sequence of reference protein
reference = "MRDTDVTVLGLGLMGQALAGAFLKDGHATTVWNRSEGKAGQLAEQGAVLASSARDAAEASPLVVVCVSDHAAVRAVLDPLGDVLAGRVLVNLTSGTSEQARATAEWAAERGITYLDGAIMAIPQVVGTADAFLLYSGPEAAYEAHEPTLRSLGAGTTYLGADHGLSSLYDVALLGIMWGTLNSFLHGAALLGTAKVEATTFAPFANRWIEAVTGFVSAYAGQVDQGAYPALDATIDTHVATVDHLIHESEAAGVNTELPRLVRTLADRALAGGQGGLGYAAMIEQFRSPS*"

# Preprocess the data
preprocessed_data = preprocessing(data_frame, 2, reference)

# Unpack the preprocessed data
single_mut_W_observed_std = preprocessed_data["Single Mutations"]["Observed std of fitness"]
sequence_double_list = preprocessed_data["Higher Order Mutations"]["Higher order mutants"]
W_observed_list = preprocessed_data["Higher Order Mutations"]["Observed fitness"]
W_observed_std_list = preprocessed_data["Higher Order Mutations"]["Observed std of fitness"]
W_expected_list = preprocessed_data["Higher Order Mutations"]["Expected fitness"]
W_expected_std_list = preprocessed_data["Higher Order Mutations"]["Expected std of fitness"]
epistatic_score_list = preprocessed_data["Higher Order Mutations"]["Epistatic score"]

candidates_tria_combined = double_mut_pos(epistatic_score_list, W_observed_list, W_expected_std_list,
                                          W_observed_std_list, sequence_double_list, reference)

# Determine all epistatic triangles for all AA positions
AA_pos_list = np.unique(candidates_tria_combined)

epistatic_triangle_list = []

for AA_pos in range(len(AA_pos_list)):
    epistatic_triangle = epistatic_interaction_double_mutation(candidates_tria_combined, AA_pos_list[AA_pos])
    if len(epistatic_triangle) > 0:
        # Append list of triangles with the new triangles one by one
        if any(isinstance(j, list) for j in epistatic_triangle):
            for triangle in range(0, len(epistatic_triangle)):
                epistatic_triangle_list.append(sorted(epistatic_triangle[triangle]))
        # Append list of triangles with the one new triangle
        else:
            epistatic_triangle_list.append(sorted(epistatic_triangle))

# List of epistatic triangles
epistatic_triangle_list = np.unique(np.array(epistatic_triangle_list), axis=0).tolist()

print(" Epistatic triangles: ", epistatic_triangle_list)

# Create epistasis double_mut_epistasis_graph given list of double mutation positions
double_mut_epistasis_graph = epistasis_graph(candidates_tria_combined)

# Plot epistasis double_mut_epistasis_graph
nx.draw(double_mut_epistasis_graph, with_labels=True, font_weight='bold')
plt.show()

# Node degree analysis (node, degree) in descending order
node_degree_list = list(map(list, sorted(double_mut_epistasis_graph.degree, key=lambda x: x[1], reverse=True)))
print(node_degree_list)

# Plot node degree list as bar chart
plt.bar(np.array(node_degree_list)[:, 0], np.array(node_degree_list)[:, 1])
plt.title("Node Degree Distribution")
plt.xlabel("Amino Acid Position")
plt.ylabel("Number of Epistatic Interactions")
plt.show()

# Max clique in double_mut_epistasis_graph (usually only a triplet)
print(nx.approximation.max_clique(double_mut_epistasis_graph))

# Plot adjacency matrix
double_mut_epistasis_graph_A = nx.adjacency_matrix(double_mut_epistasis_graph).todense()
plt.imshow(double_mut_epistasis_graph_A)
plt.show()

# Epistatic Heatmap
epistatic_heatmap = np.ones((len(reference), len(reference))) * -6

epistatic_pos_score_list = []

for double_mut in range(len(sequence_double_list)):
    epistatic_pos_score_list_element = []
    sequence = sequence_double_list[double_mut]
    _, pos, _ = call_aa_simple(reference, sequence)
    W_observed = W_observed_list[double_mut]

    # Append list element with double mutation positions and corresponding epistatic score / observed fitness
    epistatic_pos_score_list_element.append(pos)
    epistatic_pos_score_list_element.append(W_observed)

    # Append total epistatic score list with element
    epistatic_pos_score_list.append(epistatic_pos_score_list_element)

for double_mut in range(len(epistatic_pos_score_list)):
    index_pair = epistatic_pos_score_list[double_mut][0]
    score_of_index_pair = epistatic_pos_score_list[double_mut][1]
    epistatic_heatmap[index_pair[0] - 1, index_pair[1] - 1] = score_of_index_pair
    epistatic_heatmap[index_pair[1] - 1, index_pair[0] - 1] = score_of_index_pair

plt.imshow(epistatic_heatmap, cmap='viridis', interpolation='nearest')
plt.colorbar()
plt.show()

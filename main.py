import pandas as pd
from matplotlib import pyplot as plt
import networkx as nx
from analysis_utils import preprocessing, double_mut_pos, epistasis_graph, epistatic_triangles
from plotting_utils import plot_obs_fitness_heatmap, plot_node_degree_distribution

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

# Create a list of positive and combinable positions of double mutations
pos_comb_double_mut_list = double_mut_pos(epistatic_score_list, W_observed_list, W_expected_std_list,
                                          W_observed_std_list, sequence_double_list, reference)

# Determine all epistatic triangles for all AA positions
epistatic_triangle_list = epistatic_triangles(pos_comb_double_mut_list)
print(" Epistatic triangles: ", epistatic_triangle_list)

# Create epistasis double_mut_epistasis_graph given list of double mutation positions
double_mut_epistasis_graph = epistasis_graph(pos_comb_double_mut_list)

# Plot epistasis double_mut_epistasis_graph
nx.draw(double_mut_epistasis_graph, with_labels=True, font_weight='bold')
plt.show()

# Node degree analysis (node, degree) in descending order
plot_node_degree_distribution(double_mut_epistasis_graph)

# Plot adjacency matrix
double_mut_epistasis_graph_A = nx.adjacency_matrix(double_mut_epistasis_graph).todense()
plt.imshow(double_mut_epistasis_graph_A)
plt.show()

# Epistatic Heatmap
plot_obs_fitness_heatmap(reference, sequence_double_list, W_observed_list)

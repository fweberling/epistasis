import pandas as pd
from matplotlib import pyplot as plt
import networkx as nx
from analysis_utils import preprocessing, double_mut_pos, epistasis_graph, epistatic_triangles, comb_pos_mut, call_aa_simple
from plotting_utils import plot_obs_fitness_heatmap, plot_node_degree_distribution, plot_node_degree_aa_distribution, plot_mutation_distribution
import numpy as np
import itertools

# Upload input files into panda data frame
# data_frame = pd.read_csv('CPA_merge_filtered.csv')
data_frame = pd.read_csv('221111_CPA_merge_filtered.csv')
# Specify sequence of reference protein
reference = "MRDTDVTVLGLGLMGQALAGAFLKDGHATTVWNRSEGKAGQLAEQGAVLASSARDAAEASPLVVVCVSDHAAVRAVLDPLGDVLAGRVLVNLTSGTSEQARATAEWAAERGITYLDGAIMAIPQVVGTADAFLLYSGPEAAYEAHEPTLRSLGAGTTYLGADHGLSSLYDVALLGIMWGTLNSFLHGAALLGTAKVEATTFAPFANRWIEAVTGFVSAYAGQVDQGAYPALDATIDTHVATVDHLIHESEAAGVNTELPRLVRTLADRALAGGQGGLGYAAMIEQFRSPS*"

# Specify maximum order of mutations to be analysed i.e. mutations up to that order will be analysed
num_mut = 5
# Preprocess the data
preprocessed_data = preprocessing(data_frame, num_mut, reference)

############# ALL MUTATIONS
for i in range(2, num_mut + 1):
    locals()["mut_" + str(i) + "_sequence_list"] = preprocessed_data[str(i) + " Mutation"]["Sequence of mutants"]
    locals()["mut_" + str(i) + "_W_observed_list"] = preprocessed_data[str(i) + " Mutation"]["Observed fitness"]
    locals()["mut_" + str(i) + "_W_observed_std_list"] = preprocessed_data[str(i) + " Mutation"]["Observed std of " \
                                                                                                 "fitness"]
    locals()["mut_" + str(i) + "_W_expected_list"] = preprocessed_data[str(i) + " Mutation"]["Expected fitness"]
    locals()["mut_" + str(i) + "_W_expected_std_list"] = preprocessed_data[str(i) + " Mutation"]["Expected std of " \
                                                                                                 "fitness"]
    locals()["mut_" + str(i) + "_epistatic_score_list"] = preprocessed_data[str(i) + " Mutation"]["Epistatic score"]

full_mut_sequence_list = []
full_mut_W_observed_list = []
full_mut_W_observed_std_list = []
full_mut_W_expected_list = []
full_mut_W_expected_std_list = []
full_mut_epistatic_score_list = []

for mut_num_i in range(2, num_mut + 1):
    full_mut_sequence_list = full_mut_sequence_list + locals()["mut_" + str(mut_num_i) + "_sequence_list"]
    full_mut_W_observed_list = full_mut_W_observed_list + locals()["mut_" + str(mut_num_i) + "_W_observed_list"]
    full_mut_W_observed_std_list = full_mut_W_observed_std_list + locals()[
        "mut_" + str(mut_num_i) + "_W_observed_std_list"]
    full_mut_W_expected_list = full_mut_W_expected_list + locals()["mut_" + str(mut_num_i) + "_W_expected_list"]
    full_mut_W_expected_std_list = full_mut_W_expected_std_list + locals()[
        "mut_" + str(mut_num_i) + "_W_expected_std_list"]
    full_mut_epistatic_score_list = full_mut_epistatic_score_list + locals()[
        "mut_" + str(mut_num_i) + "_epistatic_score_list"]

# Plot mutation distribution for double mutations
plot_mutation_distribution(mut_2_sequence_list, reference)

# Plot mutation distribution for mutations 2 - 5
plot_mutation_distribution(full_mut_sequence_list, reference)

# Plot mutation distribution for mutations 3 - 5
mut_3_5_sequence_list = []
for mut_num_i in range(3, num_mut + 1):
    mut_3_5_sequence_list = mut_3_5_sequence_list + locals()["mut_" + str(mut_num_i) + "_sequence_list"]
plot_mutation_distribution(mut_3_5_sequence_list, reference)

# Combine two lists into combined list of double mutation positions
comb_pos_mut_pos_list, comb_pos_mut_aa_list = comb_pos_mut(full_mut_epistatic_score_list, full_mut_W_observed_list,
                                                           full_mut_W_expected_list, full_mut_W_observed_std_list,
                                                           full_mut_sequence_list, reference, 1, 1)
# Unpack list of into pairs
pos_comb_mut_edges = []
pos_comb_mut_aa = []
for higher_ord_mut in range(0, len(comb_pos_mut_pos_list)):
    higher_order_mut_list = (list(map(list, itertools.combinations(comb_pos_mut_pos_list[higher_ord_mut], 2))))
    higher_order_mut_aa_list = (list(map(list, itertools.combinations(comb_pos_mut_aa_list[higher_ord_mut], 2))))
    if len(higher_order_mut_list) == 2:
        pos_comb_mut_edges.append(higher_order_mut_list)
        pos_comb_mut_aa.append(higher_order_mut_aa_list)
    else:
        for higher_order_mut_list_ele in range(0, len(higher_order_mut_list)):
            pos_comb_mut_edges.append(higher_order_mut_list[higher_order_mut_list_ele])
            pos_comb_mut_aa.append(higher_order_mut_aa_list[higher_order_mut_list_ele])

# Create epistasis graph for all higher order mutants
higher_order_mut_epistasis_graph = epistasis_graph(pos_comb_mut_edges)

# Plot epistasis double_mut_epistasis_graph
nx.draw(higher_order_mut_epistasis_graph, with_labels=True, font_weight='bold')
plt.show()

# Node degree analysis (node, degree) in descending order
plot_node_degree_distribution(higher_order_mut_epistasis_graph)

# Node degree and amino acid distribution
pos_comb_higher_mut_pos = np.concatenate(
    (np.array(pos_comb_mut_edges)[:, 0], np.array(pos_comb_mut_edges)[:, 1].astype(int)), axis=0)
pos_comb_higher_mut_mut_aa = np.concatenate((np.array(pos_comb_mut_aa)[:, 0], np.array(pos_comb_mut_aa)[:, 1]), axis=0)

pos_comb_higher_mut_pos_aa = np.stack((pos_comb_higher_mut_pos, pos_comb_higher_mut_mut_aa), axis=1)

pos_per_aa_dict = plot_node_degree_aa_distribution(pos_comb_higher_mut_pos_aa)

############# DOUBLE MUTATIONs
# Unpack the preprocessed data
single_mut_W_observed_std = preprocessed_data["1 Mutation"]["Observed std of fitness"]
sequence_double_list = preprocessed_data["2 Mutation"]["Sequence of mutants"]
W_observed_list = preprocessed_data["2 Mutation"]["Observed fitness"]
W_observed_std_list = preprocessed_data["2 Mutation"]["Observed std of fitness"]
W_expected_list = preprocessed_data["2 Mutation"]["Expected fitness"]
W_expected_std_list = preprocessed_data["2 Mutation"]["Expected std of fitness"]
epistatic_score_list = preprocessed_data["2 Mutation"]["Epistatic score"]

# Create a list of positive and combinable positions of double mutations
pos_comb_double_mut_list_full = double_mut_pos(epistatic_score_list, W_observed_list, W_expected_std_list,
                                               W_observed_std_list, sequence_double_list, reference, 0.5, -1)
pos_comb_double_mut_list = pos_comb_double_mut_list_full[:, :2].astype(int)

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

# Node degree + amino acid distribution
pos_comb_double_mut_pos = np.concatenate(
    (pos_comb_double_mut_list_full[:, 0].astype(int), pos_comb_double_mut_list_full[:, 1].astype(int)), axis=0)
pos_comb_double_mut_aa = np.concatenate((pos_comb_double_mut_list_full[:, 2], pos_comb_double_mut_list_full[:, 3]),
                                        axis=0)
pos_comb_double_mut_pos_aa = np.stack((pos_comb_double_mut_pos, pos_comb_double_mut_aa), axis=1)

pos_per_aa_dict = plot_node_degree_aa_distribution(pos_comb_double_mut_pos_aa)

# Plot adjacency matrix
double_mut_epistasis_graph_A = nx.adjacency_matrix(double_mut_epistasis_graph).todense()
plt.imshow(double_mut_epistasis_graph_A)
plt.show()

# Fitness Heatmap
plot_obs_fitness_heatmap(reference, sequence_double_list, W_observed_list)

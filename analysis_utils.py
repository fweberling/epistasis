import numpy as np
import networkx as nx


def call_aa_simple(reference_seq, target_seq):
    """
    Takes a reference protein sequence and compares it with a possible mutant protein sequence to extract lists of
    wildtype amino acids, mutated amino acid positions, and mutated amino acid

    :param reference_seq: amino acid sequence of reference protein (string)
    :param target_seq: amino acid sequence of mutated protein (string)
    :return: mut_list_wt, mut_list_pos, mut_list_mut: List of mutated wild type amino acids, list of mutated amino acid
    positions, list of mutated amino acids
    """
    # if there are indels (not expected)

    if len(reference_seq) != len(str(target_seq)):
        print("Indel detected")
        return "indel", "indel", "indel"

    # call mutations in 3 lists if only substitutions -> dataframe would be better
    mut_list_wt = []
    mut_list_pos = []
    mut_list_mut = []
    for seq in range(0, len(reference_seq)):
        if reference_seq[seq] != target_seq[seq]:
            mut_list_wt.append(reference_seq[seq])
            mut_list_pos.append(seq + 1)
            mut_list_mut.append(target_seq[seq])

    return mut_list_wt, mut_list_pos, mut_list_mut


def preprocessing(data_frame, num_mut, reference_seq):
    """
    Pre-processing of the data to determine fitness scores, standard deviations of the fitness scores, and mutated amino
    acid positions for single mutations and the specified higher order mutation type

    :param data_frame: panda data frame containing all raw data information about the mutants
    :param num_mut: number of mutations to analyse
    :param reference_seq: string of amino acids of reference protein
    :return: preprocessed_data: nested dictionary of preprocessed data subdivided into single mutations and higher order
    mutations
    """

    # Save columns of data frame into lists
    seq_raw_list = data_frame["aa_seq"].tolist()
    seq_ids_list = data_frame.iloc[:, 0].tolist()
    nhams_list = data_frame["Nham_aa"].tolist()
    fitness_list = data_frame["fitness_comb"].tolist()
    observed_std_list = data_frame["sigma_comb"].tolist()

    # Empty lists of preprocessed data
    single_mut_positions = []
    single_mut_fitness = []
    single_mutations_mut_aa = []
    single_mutations_observed_std = []
    higher_ord_positions = []
    higher_ord_fitness = []
    higher_ord_mut_aa = []
    higher_ord_mut_observed_std = []
    pre_higher_ord_seq_list = []

    for i in range(len(seq_raw_list)):

        wt, pos, mut = call_aa_simple(reference_seq, seq_raw_list[i])
        if mut != "indel":
            n_mut = len(mut)
        else:
            n_mut = "NA"

        if n_mut == 1:
            single_mut_positions.append(pos[0])
            single_mut_fitness.append(fitness_list[i])
            single_mutations_mut_aa.append(mut[0])
            single_mutations_observed_std.append(observed_std_list[i])

        if n_mut == num_mut:
            higher_ord_positions.append(pos)
            higher_ord_fitness.append(fitness_list[i])
            higher_ord_mut_aa.append(mut)
            higher_ord_mut_observed_std.append(observed_std_list[i])
            pre_higher_ord_seq_list.append(seq_raw_list[i])

    print(f"number of variants with {num_mut} analyzed: {len(higher_ord_fitness)}")

    number_possible_epistatic_events = 0

    # Data handling epistasis (higher order mutational effects)
    W_expected_list = []
    W_observed_list = []
    W_expected_std_list = []
    W_observed_std_list = []
    epistatic_score_list = []
    W_expected_list_non_log = []
    W_observed_list_non_log = []

    # for counting epistatic hotspots
    positive_positions = []

    higher_ord_seq_list = []

    for i in range(len(higher_ord_positions)):
        # get positions of the double mutation

        positions = higher_ord_positions[i]
        mut_aas = higher_ord_mut_aa[i]
        W_observed = higher_ord_fitness[i]
        W_observed_std = higher_ord_mut_observed_std[i]

        # for every position found +1
        num_found = 0

        # list with all the single fitness_list values for this mutation
        single_fitness_list = []
        single_observed_std_list = []
        # for QC also the identity of single mutations:
        single_pos_list = []
        single_mut_aa_list = []

        for sub_pos in range(len(positions)):

            # loop through single mutant information
            for a in range(len(single_mut_positions)):
                # if match for subpositions -> append lists and numb found += 1
                if single_mut_positions[a] == positions[sub_pos] and single_mutations_mut_aa[a] == mut_aas[sub_pos]:
                    num_found += 1
                    single_fitness_list.append(single_mut_fitness[a])
                    single_observed_std_list.append(single_mutations_observed_std[a])
                    # only for control
                    single_pos_list.append(single_mut_positions[a])
                    single_mut_aa_list.append(single_mutations_mut_aa[a])

        # only if all positions found: proceed
        if num_found == len(positions):
            number_possible_epistatic_events += 1
            single_observed_std_list = np.array(single_observed_std_list)

            # product model of expected fitness:
            W_expected = sum(single_fitness_list)
            W_expected_std = np.sqrt((single_observed_std_list ** 2).sum())
            epistatic_score = W_observed - W_expected

            W_expected_list.append(W_expected)
            W_expected_std_list.append(W_expected_std)
            W_observed_list.append(W_observed)
            W_observed_std_list.append(W_observed_std)

            epistatic_score_list.append(epistatic_score)
            higher_ord_seq_list.append(pre_higher_ord_seq_list[i])
            # for model QC:
            W_observed_list_non_log.append(np.exp(W_observed))
            W_expected_list_non_log.append(np.exp(W_expected))

            """
            #additive model: 
            non_log_W_list = []
            for W in range(len(single_fitness_list)):
                non_log_W_list.append(np.exp(single_fitness_list[W]))

            W_expected = np.log(sum(non_log_W_list) - 1)
            epistatic_score = W_observed - W_expected

            W_expected_list.append(W_expected)
            W_observed_list.append(W_observed)
            epistatic_score_list.append(epistatic_score)

            """

    preprocessed_data = {
        "Single Mutations": {
            "Observed std of fitness": single_mutations_observed_std,
        },
        "Higher Order Mutations": {
            "Higher order mutants": higher_ord_seq_list,
            "Observed fitness": W_observed_list,
            "Observed std of fitness": W_observed_std_list,
            "Expected fitness": W_expected_list,
            "Expected std of fitness": W_expected_std_list,
            "Epistatic score": epistatic_score_list
        }
    }

    return preprocessed_data


def double_mut_pos(epistatic_score_list, obs_fitness_list, exp_fitness_list, obs_fitness_std_list,
                   double_mut_seq_list, reference_seq):
    """
    Creates a list of double mutation positions based on a positive epistatic effect and combinability (positive
    epistatic effect and additive single mutations)

    :param epistatic_score_list: list of expected fitness scores for each double mutant
    :param obs_fitness_list: list of observed fitness scores for each double mutant
    :param exp_fitness_list: list of standard deviations of expected fitness score for list of double mutations
    :param obs_fitness_std_list: list of standard deviations of observed fitness score for list of double mutations
    :param double_mut_seq_list: list of sequences containing double mutations (with sufficient information)
    :param reference_seq: reference sequence of protein
    :return: pos_comb_double_mut_list: list of positive and combinable double mutation positions
    """

    candidates_tria_1 = []
    candidates_tria_2 = []

    # Create lists of double mutation positions dependent on positiveness and combinability
    for seq in range(len(double_mut_seq_list)):
        # Ensure positive fitness effect and combinability
        if obs_fitness_list[seq] > obs_fitness_std_list[seq] and epistatic_score_list[seq] > - exp_fitness_list[seq]:
            sequence = double_mut_seq_list[seq]
            _, mut_pos, _ = call_aa_simple(reference_seq, sequence)

            candidates_tria_1.append(mut_pos[0])
            candidates_tria_2.append(mut_pos[1])

    # Combine two lists into combined list of double mutation positions
    candidates_tria_1_np = np.array(candidates_tria_1)
    candidates_tria_2_np = np.array(candidates_tria_2)

    pos_comb_double_mut_list = np.stack((candidates_tria_1_np, candidates_tria_2_np), axis=1).tolist()

    return pos_comb_double_mut_list


def double_mut_pos_eps(exp_fitness_list, obs_fitness_list, double_mut_seq_list, ref_seq, epistasis, add_threshold,
                       pos_threshold):
    """
    Analyses the list of sequences of double mutants for positive, negative, or no epistasis and returns a list of
    double mutation positions for each sequence

    :param exp_fitness_list: list of expected fitness scores of double mutant
    :param obs_fitness_list: list of observed fitness scores of double mutant
    :param double_mut_seq_list: list of sequences of double mutants
    :param ref_seq: reference sequence of protein
    :param epistasis: 'positive', 'negative, 'none'
    :param add_threshold: threshold of additive effect
    :param pos_threshold: threshold of positive fitness
    :return: filtered_double_mut_pos_list: list of positions of double mutations filtered according to type and
    thresholds
    """

    if epistasis != 'positive' and epistasis != 'negative' and epistasis != 'none':
        raise ValueError("Only 'positive', 'negative', or 'none' allowed as input for :param epistasis")

    if add_threshold < 0:
        raise ValueError("The value of :param add_threshold must be greater equal 0")

    if pos_threshold < 0:
        raise ValueError("The value of :param pos_threshold must be greater equal 0")

    candidates_tria_1 = []
    candidates_tria_2 = []

    # Create lists of double mutation positions dependent on epistatic or additive effect
    for k in range(len(double_mut_seq_list)):
        # Positive epistasis associated with positive outcome
        if epistasis == 'positive' and pos_threshold > 0:
            if obs_fitness_list[k] > (exp_fitness_list[k] + add_threshold) and obs_fitness_list[k] > pos_threshold:
                seq_i = double_mut_seq_list[k][0]
                _, mut_pos, _ = call_aa_simple(ref_seq, seq_i)

                candidates_tria_1.append(mut_pos[0])
                candidates_tria_2.append(mut_pos[1])
        # Positive epistasis
        elif epistasis == 'positive' and pos_threshold == 0:
            if obs_fitness_list[k] > (exp_fitness_list[k] + add_threshold):
                seq_i = double_mut_seq_list[k][0]
                _, mut_pos, _ = call_aa_simple(ref_seq, seq_i)

                candidates_tria_1.append(mut_pos[0])
                candidates_tria_2.append(mut_pos[1])
        # Negative epistasis
        elif epistasis == 'negative':
            if obs_fitness_list[k] < (exp_fitness_list[k] - add_threshold):
                seq_i = double_mut_seq_list[k][0]
                _, mut_pos, _ = call_aa_simple(ref_seq, seq_i)

                candidates_tria_1.append(mut_pos[0])
                candidates_tria_2.append(mut_pos[1])
        # No epistasis / additive effect
        elif epistasis == 'none':
            if (exp_fitness_list[k] - add_threshold) <= obs_fitness_list[k] <= (exp_fitness_list[k] + add_threshold):
                seq_i = double_mut_seq_list[k][0]
                _, mut_pos, _ = call_aa_simple(ref_seq, seq_i)

                candidates_tria_1.append(mut_pos[0])
                candidates_tria_2.append(mut_pos[1])

    # Combine two lists into combined list of double mutation positions
    candidates_tria_1_np = np.array(candidates_tria_1)
    candidates_tria_2_np = np.array(candidates_tria_2)

    filtered_double_mut_pos_list = np.stack((candidates_tria_1_np, candidates_tria_2_np), axis=1).tolist()

    return filtered_double_mut_pos_list


def epistatic_interaction_double_mutation(double_mut_pos_list, query_pos):
    """
    Analyzes a list of double mutation positions given a specific amino acid position to determine triangular
    epistatic interactions between amino acid positions

    :param double_mut_pos_list: n x 2 list of double mutation positions
    :param query_pos: amino acid position to analyze
    :return: epistatic_interaction_given_query: list of epistatic triangles given the queried amino acid position
    """

    # Convert to numpy arrayArray of combinations
    double_mut_pos_list = np.array(double_mut_pos_list)

    # Find query pairs
    idx_query_pos, _ = np.nonzero(double_mut_pos_list == query_pos)

    query_pairs = double_mut_pos_list[idx_query_pos, :]
    # print("Query pairs: ", query_pairs)

    # Find epistasis partner
    idx_1, idx_2 = np.nonzero(query_pairs != query_pos)

    epistasis_partner = query_pairs[idx_1, idx_2]
    # print("Epistasis partner: ", epistasis_partner)

    # List of epistasis partner interactions
    pairs_given_partner_list = np.zeros((1, 2))

    for partner_pos in range(len(epistasis_partner)):
        epistasis_partner_pos = epistasis_partner[partner_pos]
        idx_epistasis_partner_pos, _ = np.nonzero(double_mut_pos_list == epistasis_partner_pos)
        pairs_given_partner = double_mut_pos_list[idx_epistasis_partner_pos, :]
        pairs_given_partner_list = np.concatenate((pairs_given_partner_list, pairs_given_partner), axis=0)

    pairs_given_partner_list = pairs_given_partner_list[1:, :]
    # print("Pairs given partner: ", pairs_given_partner_list)

    # Create list of combinations of partners to compare with array
    poss_interaction_comb = np.zeros((1, 2))

    if len(epistasis_partner) == 2:
        poss_interaction_comb = np.concatenate((poss_interaction_comb, np.expand_dims(epistasis_partner, axis=0)),
                                               axis=0)
        # add other permutation
        epistasis_partner_perm = np.array([[epistasis_partner[1], epistasis_partner[0]]])
        poss_interaction_comb = np.concatenate((poss_interaction_comb, epistasis_partner_perm), axis=0)

        poss_interaction_comb = poss_interaction_comb[1:, :]

    else:

        for pos_1 in range(len(epistasis_partner)):
            for pos_2 in range(len(epistasis_partner)):
                if epistasis_partner[pos_1] != epistasis_partner[pos_2]:
                    comb = np.expand_dims(np.array([epistasis_partner[pos_1], epistasis_partner[pos_2]]), axis=0)

                    poss_interaction_comb = np.concatenate((poss_interaction_comb, comb), axis=0)

        poss_interaction_comb = poss_interaction_comb[1:, :]

    # print("Possible interaction combinations: ", poss_interaction_comb)

    # Compare

    epistatic_interaction_given_query = []

    for pair in range(len(pairs_given_partner_list)):
        for comb in range(len(poss_interaction_comb)):
            if pairs_given_partner_list[pair, 0] == poss_interaction_comb[comb, 0] and pairs_given_partner_list[
                pair, 1] == poss_interaction_comb[comb, 1]:
                epistatic_interaction_given_query.append(pairs_given_partner_list[pair, :].tolist())

    epistatic_interaction_given_query = np.array(epistatic_interaction_given_query)

    # Add query position to list of epistatic interaction if not empty
    if epistatic_interaction_given_query.size > 0:
        epistatic_interaction_given_query = np.unique(epistatic_interaction_given_query, axis=0)
        query_vector = query_pos * np.ones(len(epistatic_interaction_given_query)).reshape(-1, 1)
        epistatic_interaction_given_query = np.column_stack((query_vector, epistatic_interaction_given_query)).squeeze()
        epistatic_interaction_given_query = sorted(epistatic_interaction_given_query.tolist())

    # print("Epistatic interaction triangle: ", epistatic_interaction_given_query)
    # print(" ")

    return epistatic_interaction_given_query


def epistasis_graph(double_mut_pos_list):
    """
    Creates an epistasis double_mut_epistasis_graph given a list of positions of double mutation.

    :param double_mut_pos_list: list of lists of positions of double mutations
    :return: double_mut_epistasis_graph: epistasis double_mut_epistasis_graph
    """
    graph = nx.Graph()

    edges_list = []
    for pair in double_mut_pos_list:
        edges_list.append(tuple(pair))

    graph.add_edges_from(edges_list)

    return graph


def epistatic_triangles(higher_order_mut_positions):
    """
    Given a list of higher order mutation positions, all epistatic trianglular structures are exptracted as a list
    :param higher_order_mut_positions: list of higher order mutation positions
    :return: epistatic_triangle_list: returns a list of epistatic triangles
    """
    # Determine all epistatic triangles for all AA positions
    AA_pos_list = np.unique(higher_order_mut_positions)

    epistatic_triangle_list = []

    for AA_pos in range(len(AA_pos_list)):
        epistatic_triangle = epistatic_interaction_double_mutation(higher_order_mut_positions, AA_pos_list[AA_pos])
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

    return epistatic_triangle_list

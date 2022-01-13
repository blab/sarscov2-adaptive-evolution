from utils import get_parent
from Bio import Phylo
import numpy as np
import random



def get_total_muts_on_tree(tree):
    """
    Find the total number of mutations that occur on internal branches of the tree. Return a dictionary for nonsynonymous mutations and another for synonymous, each with gene as key and the value as the total number of mutations in that gene on internal branches of the tree
    """
    total_mutations_nonsyn = {'S1': 0, 'RdRp': 0, 'S2':0, 'N': 0, 'E':0, 'M':0, 'ORF3a':0, 'ORF6':0, 'ORF7a':0, 'ORF7b':0, 'ORF8':0, 'ORF9b':0, 'Nsp1':0, 'Nsp2':0, 'Nsp3':0, 'Nsp4':0, 'Nsp5':0, 'Nsp6':0, 'Nsp7':0, 'Nsp8':0, 'Nsp9':0, 'Nsp10':0, 'Nsp13':0, 'Nsp14':0, 'Nsp15':0, 'Nsp16':0}
    total_mutations_syn = {'S1':0, 'S2': 0, 'RdRp':0, 'E':0, 'M':0, 'N':0, 'ORF7a':0, 'Nsp6':0}

    for node in tree.find_clades(terminal=False):
        # look only at nodes with at least 3 descendents, to exclude weird stuff at tips and sequencing errors
        if len(node.get_terminals()) >=3:
            if len(node.nonsyn_at_node)!=0:
                for k,v in node.nonsyn_at_node.items():
                    total_mutations_nonsyn[k]+=v
            if len(node.syn_at_node)!=0:
                for k,v in node.syn_at_node.items():
                    total_mutations_syn[k]+=v

    return total_mutations_nonsyn, total_mutations_syn


def get_total_num_branches(tree):
    """
    Return the total number of internal branches on the get_total_muts_on_tree
    """
    total_branchs = 0
    for node in tree.find_clades(terminal=False):
        if len(node.get_terminals()) >=3:
            total_branchs+=1

    return total_branchs


def get_branch_lengths(tree, num_tips=3):
    """
    Find the length of every internal branch on the tree
    """

    root_date = tree.root.node_attrs['num_date']['value']

    # keep track of parts of paths that have already been considered
    paths_segments_already_used_all = []

    # store all wait times between branchs
    branch_lengths_all = []
    branch_names_all = []

    for node in tree.find_clades(terminal=False):

        if len(node.get_terminals()) >=num_tips:


            node_path = get_parent(tree, node)

            # find the branch length (in decimal years) between node and previous node
            path_dates = [n.node_attrs['num_date']['value'] for n in node_path]

            if len(path_dates) > 1:
                path_dates_with_root = [root_date] + path_dates
                branch_lengths = [j-i for i, j in zip(path_dates_with_root[:-1], path_dates_with_root[1:])]


                for n_index in range(len(node_path)):

                    # name of the node at the end of the branch
                    # branch length will be assigned to this node because
                    # it is the length of the branch leading to this node
                    branch_name = node_path[n_index].name
                    if n_index==0:
                        preceding_node_name = tree.root.name
                    else:
                        preceding_node_name = node_path[n_index-1].name


                    if [preceding_node_name, branch_name] not in paths_segments_already_used_all:

                        branch_length = branch_lengths[n_index]

                        if branch_length<0:
                            branch_length = 0.0

                        branch_lengths_all.append(branch_length)
                        branch_names_all.append(branch_name)
                        paths_segments_already_used_all.append([preceding_node_name, branch_name])

    return branch_names_all, branch_lengths_all


def randomize_mutations_on_tree_multinomial(tree, num_mutations, branch_lengths, branch_names):
    """
    Randomize mutation positions on the phylogeny. Do by randomly choosing X branchs to have mutations, where X is the total number of mutations of some category observed on the phylogeny (num_mutations). Branches are chosen according to a multinomial draw where the likelihood of getting a mutation is proportional to the branch length. Can result in multiple mutations on the same branch, as with empirical data
    """


    # multinomial draw to randomize mutations on tree
    branch_lengths_normalized = [x/sum(branch_lengths) for x in branch_lengths]

    num_branches = len(branch_lengths_normalized)

    hits = list(np.random.multinomial(num_mutations, branch_lengths_normalized, size=1))[0]

    muts_on_branch = dict(zip(branch_names, hits))

    for node in tree.find_clades(terminal=False):
        if len(node.get_terminals()) >=3:
        # no branch length at root
            if node.name == 'NODE_0000000' or node.name == 'NODE_0000001':
                node.random_muts = 0

            else:
                node.random_muts = muts_on_branch[node.name]
    return tree


def random_wait_times(tree):
    root_date = tree.root.node_attrs['num_date']['value']

    # keep track of parts of paths that have already been considered
    paths_segments_already_used_random = []

    # keep track of nodes (with multiple S1 mutations) where the wait times
    # of 0 (between those mutations) have already been counted
    wait_times_already_counted_at_node = []

    # store all wait times between observed mutations
    all_wait_times_random = []



    for node in tree.find_clades(terminal=False):

        if len(node.get_terminals()) >=3:



            node_path = get_parent(tree, node)

            # find the branch length (in decimal years) between node and previous node
            path_dates = [n.node_attrs['num_date']['value'] for n in node_path]

            if len(path_dates) > 0:
                path_dates_with_root = [root_date] + path_dates
                branch_lengths = [j-i for i, j in zip(path_dates_with_root[:-1], path_dates_with_root[1:])]

            random_muts_on_path = [n.random_muts for n in node_path]
            # indicies of branches that have S1 mutations
            random_mut_indicies = [i for i, num_random in enumerate(random_muts_on_path) if num_random!=0]

            # count wait time from root until first mutation(s)
            if len(random_mut_indicies) == 1:
                A_name = node_path[random_mut_indicies[0]].name
                root_name = node_path[0].name

                if [root_name, A_name] not in paths_segments_already_used_random:
                    num_muts_A = random_muts_on_path[random_mut_indicies[0]]

                    # pick date for branch A randomly between inferred node date and previous node
                    # because mutation happened somewhere in this time window, not necessarily right at the node
                    node_date_A = node_path[random_mut_indicies[0]].node_attrs['num_date']['value']

                    # if there are multiple mutations at this branch, assign each of them a different random date
                    dates_A = []
                    for mut in range(num_muts_A):
                        date_A = random.uniform(node_date_A - branch_lengths[random_mut_indicies[0]], node_date_A)
                        dates_A.append(date_A)

                    dates_A = sorted(dates_A)

                    # prepend the root date
                    dates_A = [root_date] + dates_A

                    # find wait times between this branch and the root and between mutations on this branch
                    wait_times = [j-i for i, j in zip(dates_A[:-1], dates_A[1:])]
                    all_wait_times_random+=wait_times

                    paths_segments_already_used_random.append([root_name, A_name])
                    wait_times_already_counted_at_node.append(A_name)



            # if more than 1 mutation is on this path, calculate wait time between mutations
            if len(random_mut_indicies) > 1:

                for x in range(len(random_mut_indicies) -1):

                    A_name = node_path[random_mut_indicies[x]].name
                    B_name = node_path[random_mut_indicies[x+1]].name

                    # only count wait times if this path segment hasn't already been counted
                    if [A_name, B_name] not in paths_segments_already_used_random:

                        num_muts_A = random_muts_on_path[random_mut_indicies[x]]
                        num_muts_B = random_muts_on_path[random_mut_indicies[x+1]]

                        # pick dates for branches A and B randomly between inferred node date and previous node
                        # because mutation happened somewhere in this time window, not necessarily right at the node


                        # if there are multiple mutations on these branch, assign each of them a different random date
                        node_date_A = node_path[random_mut_indicies[x]].node_attrs['num_date']['value']
                        dates_A = []
                        for mut in range(num_muts_A):
                            date_A = random.uniform(node_date_A - branch_lengths[random_mut_indicies[x]], node_date_A)
                            dates_A.append(date_A)

                        dates_A = sorted(dates_A)


                        node_date_B = node_path[random_mut_indicies[x+1]].node_attrs['num_date']['value']
                        dates_B = []
                        for mut in range(num_muts_B):
                            date_B = random.uniform(node_date_B - branch_lengths[random_mut_indicies[x+1]], node_date_B)
                            dates_B.append(date_B)

                        dates_B = sorted(dates_B)


                        # find wait times between mutations on one branch
                        # only if this node hasn't already been counted
                        if A_name not in wait_times_already_counted_at_node:
                            wait_times_on_A = [j-i for i, j in zip(dates_A[:-1], dates_A[1:])]
                            all_wait_times_random+=wait_times_on_A
                            wait_times_already_counted_at_node.append(A_name)

                        if B_name not in wait_times_already_counted_at_node:
                            wait_times_on_B = [j-i for i, j in zip(dates_B[:-1], dates_B[1:])]
                            all_wait_times_random+=wait_times_on_B
                            wait_times_already_counted_at_node.append(B_name)

                        # and between latest mutation on branch A and earliest on branch B
                        wait_time_between = dates_B[0] - dates_A[-1]
                        all_wait_times_random.append(wait_time_between)



                        # add this pair of branches to the list of path segments that have already been counted
                        paths_segments_already_used_random.append([A_name, B_name])

    return all_wait_times_random


def s1_wait_times(tree, nonsyn_syn):
    # to deal with branches with multiple mutations, say that each mutation happens at a randomly chosen time
    # between the node and previous node


    root_date = tree.root.node_attrs['num_date']['value']

    # keep track of parts of paths that have already been considered
    paths_segments_already_used_s1 = []

    # keep track of nodes (with multiple S1 mutations) where the wait times
    # of 0 (between those mutations) have already been counted
    wait_times_already_counted_at_node = []

    # store all wait times between observed mutations
    all_wait_times_s1 = []



    for node in tree.find_clades(terminal=False):

        if len(node.get_terminals()) >=3:


            node_path = get_parent(tree, node)

            # find the branch length (in decimal years) between node and previous node
            path_dates = [n.node_attrs['num_date']['value'] for n in node_path]

            if len(path_dates) > 0:
                path_dates_with_root = [root_date] + path_dates
                branch_lengths = [j-i for i, j in zip(path_dates_with_root[:-1], path_dates_with_root[1:])]

            if nonsyn_syn == 'nonsyn':
                s1_muts_on_path = [n.nonsyn_at_node['S1'] for n in node_path]
            if nonsyn_syn == 'syn':
                s1_muts_on_path = [n.syn_at_node['S1'] for n in node_path]
            # indicies of branches that have S1 mutations
            s1_mut_indicies = [i for i, num_s1 in enumerate(s1_muts_on_path) if num_s1!=0]

            # count wait time from root until first mutation(s)
            if len(s1_mut_indicies) == 1:
                A_name = node_path[s1_mut_indicies[0]].name
                root_name = node_path[0].name

                if [root_name, A_name] not in paths_segments_already_used_s1:
                    num_muts_A = s1_muts_on_path[s1_mut_indicies[0]]

                    # pick date for branch A randomly between inferred node date and previous node
                    # because mutation happened somewhere in this time window, not necessarily right at the node
                    node_date_A = node_path[s1_mut_indicies[0]].node_attrs['num_date']['value']

                    # if there are multiple mutations at this branch, assign each of them a different random date
                    dates_A = []
                    for mut in range(num_muts_A):
                        date_A = random.uniform(node_date_A - branch_lengths[s1_mut_indicies[0]], node_date_A)
                        dates_A.append(date_A)

                    dates_A = sorted(dates_A)

                    # prepend the root date
                    dates_A = [root_date] + dates_A

                    # find wait times between this branch and the root and between mutations on this branch
                    wait_times = [j-i for i, j in zip(dates_A[:-1], dates_A[1:])]
                    all_wait_times_s1+=wait_times

                    paths_segments_already_used_s1.append([root_name, A_name])
                    wait_times_already_counted_at_node.append(A_name)



            # if more than 1 mutation is on this path, calculate wait time between mutations
            if len(s1_mut_indicies) > 1:

                for x in range(len(s1_mut_indicies) -1):

                    A_name = node_path[s1_mut_indicies[x]].name
                    B_name = node_path[s1_mut_indicies[x+1]].name

                    # only count wait times if this path segment hasn't already been counted
                    if [A_name, B_name] not in paths_segments_already_used_s1:

                        num_muts_A = s1_muts_on_path[s1_mut_indicies[x]]
                        num_muts_B = s1_muts_on_path[s1_mut_indicies[x+1]]

                        # pick dates for branches A and B randomly between inferred node date and previous node
                        # because mutation happened somewhere in this time window, not necessarily right at the node


                        # if there are multiple mutations on these branch, assign each of them a different random date
                        node_date_A = node_path[s1_mut_indicies[x]].node_attrs['num_date']['value']
                        dates_A = []
                        for mut in range(num_muts_A):
                            date_A = random.uniform(node_date_A - branch_lengths[s1_mut_indicies[x]], node_date_A)
                            dates_A.append(date_A)

                        dates_A = sorted(dates_A)


                        node_date_B = node_path[s1_mut_indicies[x+1]].node_attrs['num_date']['value']
                        dates_B = []
                        for mut in range(num_muts_B):
                            date_B = random.uniform(node_date_B - branch_lengths[s1_mut_indicies[x+1]], node_date_B)
                            dates_B.append(date_B)

                        dates_B = sorted(dates_B)


                        # find wait times between mutations on one branch
                        # only if this node hasn't already been counted
                        if A_name not in wait_times_already_counted_at_node:
                            wait_times_on_A = [j-i for i, j in zip(dates_A[:-1], dates_A[1:])]
                            all_wait_times_s1+=wait_times_on_A
                            wait_times_already_counted_at_node.append(A_name)

                        if B_name not in wait_times_already_counted_at_node:
                            wait_times_on_B = [j-i for i, j in zip(dates_B[:-1], dates_B[1:])]
                            all_wait_times_s1+=wait_times_on_B
                            wait_times_already_counted_at_node.append(B_name)

                        # and between latest mutation on branch A and earliest on branch B
                        wait_time_between = dates_B[0] - dates_A[-1]
                        all_wait_times_s1.append(wait_time_between)



                        # add this pair of branches to the list of path segments that have already been counted
                        paths_segments_already_used_s1.append([A_name, B_name])



    return all_wait_times_s1


def rdrp_wait_times(tree):
    root_date = tree.root.node_attrs['num_date']['value']

    # keep track of parts of paths that have already been considered
    paths_segments_already_used_rdrp = []

    # keep track of nodes (with multiple S1 mutations) where the wait times
    # of 0 (between those mutations) have already been counted
    wait_times_already_counted_at_node_rdrp = []

    # store all wait times between observed mutations
    all_wait_times_rdrp = []


    for node in tree.find_clades(terminal=False):

        if len(node.get_terminals()) >=3:

            node_path = get_parent(tree, node)

            # find the branch length (in decimal years) between node and previous node
            path_dates = [n.node_attrs['num_date']['value'] for n in node_path]

            if len(path_dates) > 0:
                path_dates_with_root = [root_date] + path_dates
                branch_lengths = [j-i for i, j in zip(path_dates_with_root[:-1], path_dates_with_root[1:])]


            rdrp_muts_on_path = [n.nonsyn_at_node['RdRp'] for n in node_path]

            rdrp_mut_indicies = [i for i, num_rdrp in enumerate(rdrp_muts_on_path) if num_rdrp!=0]

            # count wait time from root until first mutation(s)
            if len(rdrp_mut_indicies) == 1:
                A_name = node_path[rdrp_mut_indicies[0]].name
                root_name = node_path[0].name

                if [root_name, A_name] not in paths_segments_already_used_rdrp:
                    num_muts_A = rdrp_muts_on_path[rdrp_mut_indicies[0]]

                    # pick date for branch A randomly between inferred node date and previous node
                    # because mutation happened somewhere in this time window, not necessarily right at the node
                    node_date_A = node_path[rdrp_mut_indicies[0]].node_attrs['num_date']['value']

                    # if there are multiple mutations at this branch, assign each of them a different random date
                    dates_A = []
                    for mut in range(num_muts_A):
                        date_A = random.uniform(node_date_A - branch_lengths[rdrp_mut_indicies[0]], node_date_A)
                        dates_A.append(date_A)

                    dates_A = sorted(dates_A)

                    # prepend the root date
                    dates_A = [root_date] + dates_A

                    # find wait times between this branch and the root and between mutations on this branch
                    wait_times = [j-i for i, j in zip(dates_A[:-1], dates_A[1:])]

                    all_wait_times_rdrp+=wait_times


                    paths_segments_already_used_rdrp.append([root_name, A_name])
                    wait_times_already_counted_at_node_rdrp.append(A_name)



            # only care if more than 1 mutation is on this path (to calculate wait time between mutations)
            if len(rdrp_mut_indicies) > 1:

                for x in range(len(rdrp_mut_indicies) -1):

                    A_name = node_path[rdrp_mut_indicies[x]].name
                    B_name = node_path[rdrp_mut_indicies[x+1]].name

                    # only count wait times if this path segment hasn't already been counted
                    if [A_name, B_name] not in paths_segments_already_used_rdrp:

                        num_muts_A = rdrp_muts_on_path[rdrp_mut_indicies[x]]
                        num_muts_B = rdrp_muts_on_path[rdrp_mut_indicies[x+1]]

                        # pick dates for branches A and B randomly between inferred node date and previous node
                        # because mutation happened somewhere in this time window, not necessarily right at the node


                        # if there are multiple mutations on these branch, assign each of them a different random date
                        node_date_A = node_path[rdrp_mut_indicies[x]].node_attrs['num_date']['value']
                        dates_A = []
                        for mut in range(num_muts_A):
                            date_A = random.uniform(node_date_A - branch_lengths[rdrp_mut_indicies[x]], node_date_A)
                            dates_A.append(date_A)

                        dates_A = sorted(dates_A)


                        node_date_B = node_path[rdrp_mut_indicies[x+1]].node_attrs['num_date']['value']
                        dates_B = []
                        for mut in range(num_muts_B):
                            date_B = random.uniform(node_date_B - branch_lengths[rdrp_mut_indicies[x+1]], node_date_B)
                            dates_B.append(date_B)

                        dates_B = sorted(dates_B)


                        # find wait times between mutations on one branch
                        # only if this node hasn't already been counted
                        if A_name not in wait_times_already_counted_at_node_rdrp:
                            wait_times_on_A = [j-i for i, j in zip(dates_A[:-1], dates_A[1:])]
                            all_wait_times_rdrp+=wait_times_on_A
                            wait_times_already_counted_at_node_rdrp.append(A_name)

                        if B_name not in wait_times_already_counted_at_node_rdrp:
                            wait_times_on_B = [j-i for i, j in zip(dates_B[:-1], dates_B[1:])]
                            all_wait_times_rdrp+=wait_times_on_B
                            wait_times_already_counted_at_node_rdrp.append(B_name)

                        # and between latest mutation on branch A and earliest on branch B
                        wait_time_between = dates_B[0] - dates_A[-1]
                        all_wait_times_rdrp.append(wait_time_between)

                        # add this pair of branches to the list of path segments that have already been counted
                        paths_segments_already_used_rdrp.append([A_name, B_name])


    return all_wait_times_rdrp


def remove_duplicate_mutations(mutation_list):
    """
    Function to find duplicates in mutation list
    and keep only one instance of the duplicated mutation
    """
    mutation_list = list(mutation_list)
    seen = {}
    dupes = []

    for x in mutation_list:
        if x not in seen:
            seen[x] = 1
        else:
            if seen[x] >= 1:
                dupes.append(x)
            seen[x] += 1

    nondupe_list = [x for x in seen.keys()]

    return nondupe_list, dupes


def remove_ancestral_muts(tree, list_of_muts, node):
    """
    Given a list of mutations and a node, finds whether any of these mutations occurred in ancestors of the node. Returns a list of mutations that occurred in ancestors and an edited list of mutations excluding those that occurred in ancestors
    """
    # need to check that these mutation haven't already occurred on ancestors of this branch
    node_path = get_parent(tree, node)[:-1]

    dont_use = []
    for mut in list_of_muts:
        for parent in node_path:
            if mut in parent.random_muts:
                dont_use.append(mut)


    list_of_muts_edited = [x for x in list_of_muts if x not in dont_use]

    return list_of_muts_edited, dont_use


def randomize_specific_mutations_multinomial(tree, observed_mutations, branch_lengths, branch_names):
    """
    Randomize mutation positions on the phylogeny. The position of every observed mutation is randomized according to a multinomial draw. The randomization is checked to ensure that mutations make sense phylogenetically, meaning that the same mutation doesn't occur twice on one branch or twice in the same path through the tree. Will return 'bad_randomization' if these criteria cannot be met, otherwise will return the tree with an attribute 'random_muts' listing the randomized mutations assigned to the node. Mutations are limited to internal branches of the tree with at least 15 descending tips (because this is where the observed_mutations occured).

    Each observed mutation will occur the same number of times on the randomized tree as it does on the empirical tree.
    """

    # normalize the branch length to determine likelihood of mutation occurring on this branch
    branch_lengths_normalized = [x/sum(branch_lengths) for x in branch_lengths]

    num_branches = len(branch_lengths_normalized)

    # randomly throw mutations on branches according to multinomial draw
    # this decides how many mutations is on each branch (including 0)
    hits = list(np.random.multinomial(len(observed_mutations), branch_lengths_normalized, size=1))[0]

    # dictionary of how many mutations should be on each branch, by branch name
    muts_on_branch = dict(zip(branch_names, hits))

    # shuffle the bag of mutations
    shuffled_bag_of_muts = np.random.choice(observed_mutations, len(observed_mutations), replace = False)

    for node in tree.find_clades(terminal=False):
        node.random_muts = []
        if len(node.get_terminals()) >=15:
        # no branch length at root
            if node.name == 'NODE_0000000' or node.name == 'NODE_0000001':
                node.random_muts = []
            else:
                # how many muts should be on this branch
                num_muts_on_this_branch = muts_on_branch[node.name]

                if num_muts_on_this_branch == 0:
                    node.random_muts = []

                elif num_muts_on_this_branch != 0:
                    bag_before = len(shuffled_bag_of_muts)

                    # randomly choose that many mutations to throw on this node. Without replacement
                    randomized_muts = shuffled_bag_of_muts[-num_muts_on_this_branch:]
                    # then remove the chosen muts from the bag
                    shuffled_bag_of_muts = shuffled_bag_of_muts[:-num_muts_on_this_branch]


                    add_back_to_bag = []

                    dont_use = ['filler']
                    dupes = ['filler']

                    while (len(dont_use)+len(dupes))!=0:
                        # make sure these mutations don't occur in ancestors of this node
                        randomized_muts, dont_use = remove_ancestral_muts(tree, randomized_muts, node)
                        if len(dont_use) != 0:
                            randomized_muts+=(x for x in shuffled_bag_of_muts[-len(dont_use):])
                            # remove the newly chosen muts from the bag
                            shuffled_bag_of_muts = shuffled_bag_of_muts[:-len(dont_use)]

                            # add unused mutations back to bag
                            add_back_to_bag+=dont_use

                        randomized_muts, dupes = remove_duplicate_mutations(randomized_muts)
                        if len(dupes) != 0:
                            randomized_muts+=(x for x in shuffled_bag_of_muts[-len(dupes):])
                            # remove the newly chosen muts from the bag
                            shuffled_bag_of_muts = shuffled_bag_of_muts[:-len(dupes)]
                            # add unused mutations back to bag
                            add_back_to_bag+=dupes

                    # add the duplicates back into the bag
                    shuffled_bag_of_muts = np.append(shuffled_bag_of_muts, np.array(add_back_to_bag))

                    # and shuffle the bag
                    random.shuffle(shuffled_bag_of_muts)

                    # add the random mutations to this node
                    node.random_muts = randomized_muts


    if len(shuffled_bag_of_muts) == 0:
        return tree
    else:
        return 'bad_randomization'


def randomize_one_mut_on_tree_multinomial(tree, num_occurrences):
    """
    Given the number of times a mutation occurs on the phylogeny,
    randomly throw it on the phylogeny using a multinomial draw.
    Return a list of branches with this mutation or return
    "bad_randomization" if the mutation occurs twice on the same path
    """

    branch_names, branch_lengths = get_branch_lengths(tree, num_tips=15)

    # multinomial draw to randomize mutations on tree
    branch_lengths_normalized = [x/sum(branch_lengths) for x in branch_lengths]

    hits = list(np.random.multinomial(num_occurrences, branch_lengths_normalized, size=1))[0]

    branchs_with_mut = dict(zip(branch_names, hits))
    branchs_with_mut_list = [k for k,v in branchs_with_mut.items() if v!=0]

    parents_of_all_hits = []

    for node in tree.find_clades():
        if node.name in branchs_with_mut_list:
            parents = get_parent(tree, node)[:-1]
            parent_names= [p.name for p in parents]
            parents_of_all_hits+=parent_names

    for node in tree.find_clades():
        if node.name in branchs_with_mut_list:
            if node.name in parents_of_all_hits:
                branchs_with_mut_list = 'bad_randomization'

    return branchs_with_mut_list

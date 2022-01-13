import argparse
import json
import requests
from utils import sars2_genome_info, get_parent, add_syn_mut_attribute, add_mut_at_node_attr
from utils_randomization import get_total_muts_on_tree, get_branch_lengths
from augur.utils import json_to_tree

def get_tree_muts_and_branch_lengths(tree_url):
    """
    Get information about number of mutations and branch lengths needed for randomizations
    """

    # download tree
    tree_json = requests.get(tree_url).json()

    #Put tree in Bio.Phylo format
    tree = json_to_tree(tree_json)

    # add synonymous mutations as attribute of tree
    tree = add_syn_mut_attribute(tree)

    # add attribute listing number of mutations that occurred in different genes at that node
    tree = add_mut_at_node_attr(tree)

    # count all nonsyn and syn mutations
    total_mutations_nonsyn, total_mutations_syn = get_total_muts_on_tree(tree)

    # find all branch lengths and names of corresponding branchs
    branch_names_all, branch_lengths_all = get_branch_lengths(tree)

    return total_mutations_nonsyn, total_mutations_syn, branch_names_all, branch_lengths_all


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree-url", required=True, help="url link to SARS-CoV-2 tree")
    parser.add_argument("--mut-counts-output", required=True, help="JSON with counts of mutations to all genes")
    parser.add_argument("--branches-output", required=True, help="JSON with all branch lengths and names")

    args = parser.parse_args()

    total_mutations_nonsyn, total_mutations_syn, branch_names_all, branch_lengths_all = get_tree_muts_and_branch_lengths(args.tree_url)

    # store mutation information in a json
    mutation_counts = {'nonsyn':total_mutations_nonsyn, 'syn': total_mutations_syn}

    with open(args.mut_counts_output, 'w') as outfile:
        json.dump(mutation_counts, outfile)

    # store branch information in a json
    branch_info = {'names':branch_names_all, 'lengths': branch_lengths_all}

    with open(args.branches_output, 'w') as outfile:
        json.dump(branch_info, outfile)

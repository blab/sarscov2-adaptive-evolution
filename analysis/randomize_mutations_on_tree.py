import argparse
import requests
import json
from augur.utils import json_to_tree
from utils import sars2_genome_info, get_parent
from utils_randomization import randomize_mutations_on_tree_multinomial
from scipy import stats
import pandas as pd

# read in information about gene lengths
reference_gene_locations, reference_gene_codon, gene_lengths_aa = sars2_genome_info()

def run_growth_randomizations(tree_url, total_mutations_nonsyn, total_mutations_syn, branch_lengths_all, branch_names_all, gene, iteration):
    """
    Randomizes all observed mutations across the tree and calculates the r-value for the correlation between growth rate and mutation accumulation
    """
    #Download tree json
    tree_json = requests.get(tree_url).json()

    #Put tree in Bio.Phylo format
    tree = json_to_tree(tree_json)

    # check if mutations are synonymous
    if 'syn' in gene:
        gene, nonsyn_syn = gene.split('_')

        tree = randomize_mutations_on_tree_multinomial(tree, total_mutations_syn[gene],
                                                       branch_lengths_all, branch_names_all)

    else:
        nonsyn_syn = 'nonsyn'
        tree = randomize_mutations_on_tree_multinomial(tree, total_mutations_nonsyn[gene],
                                                       branch_lengths_all, branch_names_all)


    # keep track of logistic growth rate for this iteration
    muts_information = []
    # count accumulation of mutations on path
    for node in tree.find_clades(terminal=False):

        # only care if it has a logistic growth rate
        if "logistic_growth" in node.node_attrs:
            logistic_growth = node.node_attrs["logistic_growth"]["value"]

            #Find all parents of the node
            parents = get_parent(tree, node)

            #Find mutations that occur in the parents, and at node
            parents_random_muts = 0

            for parent in parents:
                if hasattr(parent, "random_muts"):
                    parents_random_muts+=parent.random_muts

            muts_information.append({'muts_per_codon': parents_random_muts/gene_lengths_aa[gene],
                                     'logistic_growth': logistic_growth})




    muts_information_df = pd.DataFrame(muts_information)
    slope, intercept, r_value, p_value, std_err = stats.linregress(muts_information_df['logistic_growth'], muts_information_df['muts_per_codon'])

    growth_stats_dict = {'data': 'randomized', 'iteration':iteration, 'gene': gene, 'nonsyn_syn': nonsyn_syn, 'r_value':r_value}

    return growth_stats_dict



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree-url", required=True, help="url for SARS-CoV-2 tree")
    parser.add_argument("--mutation-counts", required=True, help="total number of mutations on internal branches of tree")
    parser.add_argument("--branch-info", required=True, help="names and lengths of internal branches")
    parser.add_argument("--gene", required=True, help="gene to test")
    parser.add_argument("--iteration", required=True, help="iteration of randomization")
    parser.add_argument("--output", required=True, help="JSON with randomization results for a single iteration")

    args = parser.parse_args()

    # get counts of mutations on internal branches
    with open(args.mutation_counts) as muts_json:
        muts_data = json.load(muts_json)

    # get lengths and names of internal branches
    with open(args.branch_info) as branches_json:
        branches_data = json.load(branches_json)

    growth_stats_dict = run_growth_randomizations(args.tree_url, muts_data['nonsyn'], muts_data['syn'], branches_data['lengths'], branches_data['names'], args.gene, args.iteration)

    with open(args.output, "w") as output_handle:
        json.dump(growth_stats_dict, output_handle)

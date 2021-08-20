from utils import sars2_genome_info, get_parent
from Bio import Phylo
from collections import Counter
import pandas as pd

# store gene lengths dictionary
reference_gene_locations, reference_gene_codon, gene_lengths_aa = sars2_genome_info()

def count_independent_occurrences(tree, mutation_type='aa'):

    """
    Function to determine all mutations that occur on internal branches of the phylogeny.
    Only branches with at least 15 descendents are used in order to reduce the impact of
    sequencing error and phylogenetic reconstruction (seen mostly at tips) on this analysis

    `mutation_type` can be set to 'aa' (to look at specific amino acid substitution)
    or 'site' (to look at any nonsynonymous substitution that happens at a specific residue)

    returns `all_mutations_unique`: a list of all mutations seen on internal branches and
    `independent_occurrences_all_mutations`: a dictionary with each observed mutation as a key,
    and the number of times that mutation occurs on internal branches as the value of that key.
    """

    all_mutations = []


    #only look at mutations on internal branches
    for node in tree.find_clades(terminal=False):
        # only consider mutations on branches that give rise to a clade of at least 15 tips
        if len(node.get_terminals())>=15:

            if hasattr(node, 'branch_attrs'):
                for gene, mut_list in node.branch_attrs["mutations"].items():
                    # only interested in convergent evolution of nonsyn muts for now
                    if gene!= 'nuc':
                        for mut in mut_list:

                            if mutation_type == 'aa':
                                # look at specific nonsyn muts
                                gene_mutation = f'{gene}:{mut}'


                            elif mutation_type == 'site':
                                # look at mutation to a codon (without regard to identity of mutation)
                                gene_mutation = f'{gene}:{mut[1:-1]}'


                            # append mutation to a list of all observed mutations
                            # exclude if mutation is in stop codon position of gene
                            if int(mut[1:-1])!=gene_lengths_aa[gene]:
                                all_mutations.append(gene_mutation)


    # Find all mutations that occurred on internal branches of the phylogeny
    all_mutations_unique = list(set(all_mutations))

    # Count how many times each mutation occurred on internal branches of the phylogeny
    independent_occurrences_all_mutations = Counter(all_mutations)


    return all_mutations_unique, independent_occurrences_all_mutations


def logistic_growth_of_clades_w_mut(tree, all_mutations_unique, mutation_type='aa'):

    """
    Find logistic growth rate of each clades where mutation occurs.

    Returns a dictionary with mutation as key and the value is a list of the logistic growth rates at each occurrence of the mutation
    """

    #initialize dict for storing all growth rates observed in clades with {key} mutation
    growth_dict = {k:[] for k in all_mutations_unique}


    for node in tree.find_clades(terminal=False):

        # only considering branches that give rise to a clade of at least 15
        if len(node.get_terminals())>=15:

            logistic_growth = None
            if "logistic_growth" in node.node_attrs:
                logistic_growth = node.node_attrs["logistic_growth"]["value"]

            if hasattr(node, 'branch_attrs'):
                for gene, mut_list in node.branch_attrs["mutations"].items():
                    if gene!= 'nuc':
                        for mut in mut_list:

                            if mutation_type == 'aa':
                                # look at specific nonsyn muts
                                gene_mutation = f'{gene}:{mut}'


                            elif mutation_type == 'site':
                                # look at mutation to a codon (without regard to identity of mutation)
                                gene_mutation = f'{gene}:{mut[1:-1]}'

                            # exclude mutations to last codon of gene (stop codon)
                            if int(mut[1:-1])!=gene_lengths_aa[gene]:

                                if logistic_growth!=None:
                                    growth_dict[gene_mutation].append(logistic_growth)


    # get rid of empty entries
    growth_dict = {k:v for k,v in growth_dict.items() if v}


    return growth_dict


def calc_mean_rate(growth_dict):
    """
    Calculate the mean growth rate of all clades that have mutation
    """
    growth_by_precedingmut_mean = {k:(sum(v)/len(v)) for k,v in growth_dict.items()}
    return growth_by_precedingmut_mean


def convergent_evo_dataframe(tree, all_mutations_unique, independent_occurrences_all_mutations, mutation_type='aa'):
    """
    Make a dataframe where each row is a mutations and columns contain number of independent occurrences, a list of growth rates for all occurrences, and the average growth rate
    """

    growth_dict = logistic_growth_of_clades_w_mut(tree, all_mutations_unique)
    growth_by_precedingmut_mean = calc_mean_rate(growth_dict)


    convergent_evo_list = []

    for mut, occurrences in independent_occurrences_all_mutations.items():

        if mut in growth_by_precedingmut_mean.keys():
            avg_growth = growth_by_precedingmut_mean[mut]
        else:
            avg_growth = None

        if mut in growth_dict.keys():
            rate_list = growth_dict[mut]
        else:
            rate_list = []



        # for deletions that span multiple consecutive residues, condense this into one row
        if mut == 'ORF1a:S3675-':
            mut = 'ORF1a:3675-3677del'


        elif mut == 'S:H69-':
            mut = 'S:69/70del'

        if mut not in ['ORF1a:G3676-', 'ORF1a:F3677-', 'S:V70-']:
#             for o in rate_list:
            convergent_evo_list.append({'mutation': mut, 'independent_occurrences': occurrences,
                                        'avg_growth': avg_growth, 'rate_list': rate_list
                                       })

    convergent_evo_df = pd.DataFrame(convergent_evo_list)

    return convergent_evo_df


def mean_rate_on_tree(tree):
    """
    Finds the mean growth rate of all clades on the tree.
    Limit nesting by removing duplicate growth rate values
    """

    all_growth_rates = []
    for node in tree.find_clades(terminal=False):

        # still only considering branches that give rise to a clade of at least 15
        if len(node.get_terminals())>=15:
            if "logistic_growth" in node.node_attrs:
                all_growth_rates.append(node.node_attrs["logistic_growth"]["value"])

    all_growth_rates_limit_nesting = list(set(all_growth_rates))
    mean_growth_rate = sum(all_growth_rates_limit_nesting)/len(all_growth_rates_limit_nesting)

    return mean_growth_rate


def rates_randomized_mut(all_mutations_unique, randomized_tree):

    """
    Find logistic growth rate of each clades where randomized mutation occurs
    """

    #initialize dict for storing all growth rates observed in clades with {key} mutation
    growth_dict = {k:[] for k in all_mutations_unique}

    for node in randomized_tree.find_clades():

        # still only considering branches that give rise to a clade of at least 15
        if len(node.get_terminals())>=15:



            logistic_growth = None
            if "logistic_growth" in node.node_attrs:
                logistic_growth = node.node_attrs["logistic_growth"]["value"]


            for gene_mutation in node.random_muts:

                if logistic_growth!=None:
                    growth_dict[gene_mutation].append(logistic_growth)

    # get rid of empty entries
    growth_dict = {k:v for k,v in growth_dict.items() if v}


    return growth_dict

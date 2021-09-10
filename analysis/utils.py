import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import Phylo


"""
Read in the reference file and generate dictionaries that store sequence and position information for the functions below
"""
# make dictionary with gene name as key and reference sequence of that gene as value
reference_sequence_aa = {}
reference_sequence_nt = {}

# make dictionary giving gene by genomic location
reference_gene_locations = {}

# make dictionary saying what codon within the gene a certain genomic location falls within
# and whether the mutation is at pos 0, 1 or 2 within codon
reference_gene_codon = {}

for record in SeqIO.parse(open("reference_seq_edited.gb","r"), "genbank"):
    genome_seq = record.seq
    for feature in record.features:
        if feature.type == 'CDS':
            # allow RdRp to overwrite Orf1a and Orf1b,
            # to take care of changed reading frame due to  ribosome slippage
            # S1 and S2 will also overwrite spike
            for pos in range(int(feature.location.start), int(feature.location.end)):
                reference_gene_locations[pos] = feature.qualifiers['gene'][0]
                codon_num = math.floor((pos-feature.location.start)/3)
                pos_in_codon = ((pos-feature.location.start)-codon_num*3)
                reference_gene_codon[pos] = (codon_num, pos_in_codon)

            gene_seq = feature.location.extract(record.seq)
            reference_sequence_nt[feature.qualifiers['gene'][0]] = gene_seq
            gene_seq_aa = gene_seq.translate()
            reference_sequence_aa[feature.qualifiers['gene'][0]] = gene_seq_aa

# make dictionary with length of each gene
gene_lengths_aa = {k:len(v) for k,v in reference_sequence_aa.items()}
# Nsps are part of ORF1a, but not labeled in the reference file
# manually add their lengths here (from https://www.ncbi.nlm.nih.gov/gene/43740578)
gene_lengths_aa['Nsp1'] = 180
gene_lengths_aa['Nsp2'] = 638
gene_lengths_aa['Nsp3'] = 1945
gene_lengths_aa['Nsp4'] = 500
gene_lengths_aa['Nsp5'] = 306
gene_lengths_aa['Nsp6'] = 290
gene_lengths_aa['Nsp7'] = 83
gene_lengths_aa['Nsp8'] = 198
gene_lengths_aa['Nsp9'] = 113
gene_lengths_aa['Nsp10'] = 139
gene_lengths_aa['Nsp13'] = 601
gene_lengths_aa['Nsp14'] = 527
gene_lengths_aa['Nsp15'] = 346
gene_lengths_aa['Nsp16'] = 298


def sars2_genome_info():
    """
    Returns information about the genome.

    'reference_gene_locations' is a dictionary with nucleotide position as key and gene it belongs to as value
    'reference_gene_codon' is a dictionary with nucleotide position as key and the codon it belongs to as a value. The codon is listed as a tuple in the format (codon number relative to start of gene, position within codon). Both numbers are 0-parent_diff_zero_based
    'gene_lengths_aa' is a dictionary with gene as key and length of that gene (in amino acids) as the value
    """
    return reference_gene_locations, reference_gene_codon, gene_lengths_aa


def nuc_changes_from_ref(muts_on_path):
    """
    From all the of the nucleotide changes that have occurred on the path from root to branch,
    find the most recent nuc mutation at each site (giving the genotype at the branch)
    """

    final_muts_from_ref = {}

    # overwrites genotypes at pos in historical order
    for x in muts_on_path:
        x_pos = int(x[1:-1])
        final_muts_from_ref[x_pos] = x[-1]


    return final_muts_from_ref


def determine_synonymous(nuc_muts_on_branch, parent_diffs_from_ref):
    """
    Check every nucleotide mutation that occurred on a branch to determine whether or not it is synonymous.

    For each node, all nucleotide mutations that occurred in parents of the node are applied to the reference sequence to give the genome prior to this node. Then, each nucleotide mutation at the node is made to the appropriate codon from this genome and determined to be synonymous or nonsynonymous.

    Returns a dictionary of synonymous mutations where the key is a gene and the value is a list of synonymous mutations in this gene.
    """
    parent_diffs_pos = [int(k) for k,v in parent_diffs_from_ref.items()]


    # make dictionary of synonymous (and noncoding) mutations to add to tree
    syn_muts = {}

    # don't care about deletions because they are obviously not synonymous
    for mut in nuc_muts_on_branch:
        if mut[-1]!= '-' and mut[0]!='-':
            mut_pos = int(mut[1:-1])
            # find what gene this mut happens in
            if (mut_pos-1) in reference_gene_locations.keys():
                mut_gene = reference_gene_locations[mut_pos-1]
                mut_codon_num = reference_gene_codon[mut_pos-1][0]
                mut_codon_pos = reference_gene_codon[mut_pos-1][1]

                # find the reference sequence of the codon this mutation occurs in
                codon_ref_aa = reference_sequence_aa[mut_gene][mut_codon_num]

                codon_ref_nt = reference_sequence_nt[mut_gene][(mut_codon_num*3):(mut_codon_num*3+3)]

                # check if a mutation occurred within the same codon in a parent
                # and if so, change the reference codon sequence accordingly,
                # to tell whether the mutation at this branch is synonymous or not
                codon_genome_pos = list(range((mut_pos-1-mut_codon_pos),(mut_pos-1-mut_codon_pos+3)))

                parent_codon = codon_ref_nt
                for parent_diff in parent_diffs_pos:
                    parent_diff_zero_based = parent_diff-1
                    if parent_diff_zero_based in codon_genome_pos:
                        parent_diff_pos = codon_genome_pos.index(parent_diff_zero_based)
                        parent_codon = MutableSeq(str(codon_ref_nt))
                        parent_codon[parent_diff_pos] = parent_diffs_from_ref[parent_diff]
                        parent_codon = parent_codon.toseq()


                codon_mutated = MutableSeq(str(parent_codon))
                #if deletion (or seq error) has happened at neighboring nucleotide
                if '-' in codon_mutated:
                    pass
                else:
                    codon_mutated[mut_codon_pos] = mut[-1]
                    codon_mutated = codon_mutated.toseq()
                    codon_mutated_translation = codon_mutated.translate()

                    if str(codon_ref_aa) == str(codon_mutated_translation):
                        if mut_gene in syn_muts.keys():
                            syn_muts[mut_gene] += [mut]
                        else:
                            syn_muts[mut_gene] = [mut]



            else:
                if 'noncoding' in syn_muts.keys():
                    syn_muts['noncoding'] += [mut]
                else:
                    syn_muts['noncoding'] = [mut]

    return syn_muts


def get_parent(tree, child_clade):
    """
    Function that returns the path from root to specified clade
    """
    node_path = tree.get_path(child_clade)
    return node_path


def add_syn_mut_attribute(tree):
    """
    For each node on the tree, add a node attribute 'syn_muts', which is a dictionary with genes as keys and the value is a list of all synonymous mutations that occur at that node within the gene
    """

    for node in tree.find_clades():

        node.node_attrs['syn_muts'] = {}

        # only care if this branch has some nucleotide mutations
        if hasattr(node, 'branch_attrs'):
            if 'nuc' in node.branch_attrs['mutations']:

                nuc_muts_on_branch = node.branch_attrs['mutations']['nuc']

                node_path = get_parent(tree, node)

                nucleotide_mut_path = []

                # find all nucleotide mutations that happened in parents,
                # in case they affect codons mutated on this branch
                for parent in node_path[-1]:
                    if hasattr(parent, 'branch_attrs'):
                        if 'nuc' in parent.branch_attrs['mutations']:
                            nucleotide_mut_path+=parent.branch_attrs['mutations']['nuc']

                parent_diffs_from_ref = nuc_changes_from_ref(nucleotide_mut_path)

                syn_muts_dict = determine_synonymous(nuc_muts_on_branch, parent_diffs_from_ref)

                node.node_attrs['syn_muts'] = syn_muts_dict
    return tree


def consolidate_deletions(mutation_list):
    """
    For deletion mutations, consider adjacent sites as part of the same deletion. Return a count of all mutations, where each consolidated deletion is counted as one mutation
    """
    mutation_count = len([x for x in mutation_list if x[-1]!='-'])
    deletions_only = sorted([int(x[1:-1]) for x in mutation_list if x[-1]=='-'])

    # if there are deletions, count a run of consecutive sites as a single deletion/mutation
    if len(deletions_only) != 0:
        mutation_count+=1

        deletion_tracker = deletions_only[0]
        for deleted_pos in deletions_only[1:]:
            if deleted_pos == deletion_tracker+1:
                pass
            else:
                mutation_count+=1
            deletion_tracker = deleted_pos
    return mutation_count


def consolidate_deletions_2(mutation_list):
    """
    Consolidates deletions at adjacent sites into one deletion event. Returns a list of all mutations, including with consolidated deletions
    """

    without_deletions = [x for x in mutation_list if x[-1]!='-' and x[0]!='-']
    #consolidate deletions and reversions
    deletions_only = [x for x in mutation_list if x[-1]=='-' or x[0]=='-']
    deletions_only.sort(key=lambda x:x[1:-1])


    #keep track of start of separate deletions
    separate_deletions = []

    # if there are deletions, count a run of consecutive sites as a single deletion/mutation
    if len(deletions_only) != 0:
        separate_deletions.append(deletions_only[0])

        deletion_tracker = int(deletions_only[0][1:-1])

        for deletion in deletions_only[1:]:

            deleted_pos = int(deletion[1:-1])
            if deleted_pos == deletion_tracker+1:
                pass
            else:
                separate_deletions.append(deletion)
            deletion_tracker = deleted_pos

    consolidated_mutation_list = separate_deletions + without_deletions

    return consolidated_mutation_list


def remove_reversions(mutation_list):
    """
    If site mutates and then reverts, do not count this in the mutation tally.
    If site mutates and then mutates again (but not a reversion), count only the second mutation
    """
    mutation_list_pos = [int(x[1:-1]) for x in mutation_list]
    sites_mutated_twice = set([x for x in mutation_list_pos if mutation_list_pos.count(x) > 1])
    # find if twice-mutated site was a reversion
    for site in sites_mutated_twice:
        muts_at_site = [mut for mut in mutation_list if int(mut[1:-1])==site]

        # if site reverts, remove all mutations at this site
        if muts_at_site[0][0] == muts_at_site[-1][-1]:
            for mut in range(len(muts_at_site)):
                mutation_list.remove(muts_at_site[mut])
        # if the site mutates multiple times, but doesn't revert, keep last mutation
        else:
            for mut in range(len(muts_at_site)-1):
                mutation_list.remove(muts_at_site[mut])
    return mutation_list


def add_mut_accumulation_attr(tree):
    """
    For each node on the tree, count the number of mutations that have happened along the path from root to node (including this node). Separate the counts by gene and store them as an attribute of the node. Deletions are counted with nonsynonymous mutations. A deletion of several adjacent amino acids is counted as one mutation event
    """

    for node in tree.find_clades():

        #Find all parents of the node
        parents = get_parent(tree, node)

        #Find mutations that occur in the parents
        parents_spike_muts = []
        parents_s1_muts = []
        parents_s2_muts = []
        parents_rdrp_muts = []

        parents_n_muts = []
        parents_e_muts = []
        parents_m_muts = []
        parents_orf3a_muts = []
        parents_orf6_muts = []
        parents_orf7a_muts = []
        parents_orf7b_muts = []
        parents_orf8_muts = []
        parents_orf9b_muts = []

        parents_nsp1_muts = []
        parents_nsp2_muts = []
        parents_nsp3_muts = []
        parents_nsp4_muts = []
        parents_nsp5_muts = []
        parents_nsp6_muts = []
        parents_nsp7_muts = []
        parents_nsp8_muts = []
        parents_nsp9_muts = []
        parents_nsp10_muts = []
        parents_nsp13_muts = []
        parents_nsp14_muts = []
        parents_nsp15_muts = []
        parents_nsp16_muts = []

        parents_s1_syn = []

        for parent in parents:
            if hasattr(parent, "branch_attrs") and "mutations" in parent.branch_attrs:
                if "S" in parent.branch_attrs["mutations"]:
                    parents_spike_muts+=parent.branch_attrs["mutations"]["S"]
                    for mut in parent.branch_attrs["mutations"]["S"]:
                        # nextstrain calls pos 13 in sigpep, not S1
                        if int(mut[1:-1]) in range(14,686):
                            parents_s1_muts+=[mut]
                        elif int(mut[1:-1]) in range(687,1274):
                            parents_s2_muts+=[mut]
                #find RdRp muts
                #and Nsp4 and 6 muts
                if "ORF1a" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["ORF1a"]:
                        if int(mut[1:-1]) in range(4393,4401):
                            #renumber mut according to rdrp protein
                            rdrp_mut = f'{mut[0]}{int(mut[1:-1])-4393}{mut[-1]}'
                            parents_rdrp_muts+=[rdrp_mut]

                        elif int(mut[1:-1]) in range(1,180):
                            nsp1_mut = f'{mut[0]}{int(mut[1:-1])-1}{mut[-1]}'
                            parents_nsp1_muts+=[nsp1_mut]
                        elif int(mut[1:-1]) in range(181,818):
                            nsp2_mut = f'{mut[0]}{int(mut[1:-1])-181}{mut[-1]}'
                            parents_nsp2_muts+=[nsp2_mut]
                        elif int(mut[1:-1]) in range(819,2763):
                            nsp3_mut = f'{mut[0]}{int(mut[1:-1])-819}{mut[-1]}'
                            parents_nsp3_muts+=[nsp3_mut]
                        elif int(mut[1:-1]) in range(2764,3263):
                            nsp4_mut = f'{mut[0]}{int(mut[1:-1])-2764}{mut[-1]}'
                            parents_nsp4_muts+=[nsp4_mut]
                        elif int(mut[1:-1]) in range(3264,3569):
                            nsp5_mut = f'{mut[0]}{int(mut[1:-1])-3264}{mut[-1]}'
                            parents_nsp5_muts+=[nsp5_mut]
                        elif int(mut[1:-1]) in range(3570,3859):
                            # exclude this ancestral mut
                            if mut!= 'K3833N':
                                #renumber mut according to nsp6 protein
                                nsp6_mut = f'{mut[0]}{int(mut[1:-1])-3570}{mut[-1]}'
                                parents_nsp6_muts+=[nsp6_mut]
                        elif int(mut[1:-1]) in range(3860,3942):
                            nsp7_mut = f'{mut[0]}{int(mut[1:-1])-3860}{mut[-1]}'
                            parents_nsp7_muts+=[nsp7_mut]
                        elif int(mut[1:-1]) in range(3943,4140):
                            nsp8_mut = f'{mut[0]}{int(mut[1:-1])-3943}{mut[-1]}'
                            parents_nsp8_muts+=[nsp8_mut]
                        elif int(mut[1:-1]) in range(4141,4253):
                            nsp9_mut = f'{mut[0]}{int(mut[1:-1])-4141}{mut[-1]}'
                            parents_nsp9_muts+=[nsp9_mut]
                        elif int(mut[1:-1]) in range(4254,4392):
                            nsp10_mut = f'{mut[0]}{int(mut[1:-1])-4254}{mut[-1]}'
                            parents_nsp10_muts+=[nsp10_mut]
                if "ORF1b" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["ORF1b"]:
                        if int(mut[1:-1]) in range(1,923):
                            #renumber mut according to rdrp protein
                            rdrp_mut = f'{mut[0]}{int(mut[1:-1])+9}{mut[-1]}'
                            parents_rdrp_muts+=[rdrp_mut]
                        elif int(mut[1:-1]) in range(924,1524):
                            nsp13_mut = f'{mut[0]}{int(mut[1:-1])-924}{mut[-1]}'
                            parents_nsp13_muts+=[nsp13_mut]
                        elif int(mut[1:-1]) in range(1525,2051):
                            nsp14_mut = f'{mut[0]}{int(mut[1:-1])-1525}{mut[-1]}'
                            parents_nsp14_muts+=[nsp14_mut]
                        elif int(mut[1:-1]) in range(2052,2397):
                            nsp15_mut = f'{mut[0]}{int(mut[1:-1])-2052}{mut[-1]}'
                            parents_nsp15_muts+=[nsp15_mut]
                        elif int(mut[1:-1]) in range(2398,2695):
                            nsp16_mut = f'{mut[0]}{int(mut[1:-1])-2398}{mut[-1]}'
                            parents_nsp16_muts+=[nsp16_mut]
                # find N muts
                if "N" in parent.branch_attrs["mutations"]:
                    parents_n_muts+=parent.branch_attrs["mutations"]["N"]
                # find E muts
                if "E" in parent.branch_attrs["mutations"]:
                    parents_e_muts+=parent.branch_attrs["mutations"]["E"]
                # find M muts
                if "M" in parent.branch_attrs["mutations"]:
                    parents_m_muts+=parent.branch_attrs["mutations"]["M"]
                # find ORF3a muts
                if "ORF3a" in parent.branch_attrs["mutations"]:
                    parents_orf3a_muts+=parent.branch_attrs["mutations"]["ORF3a"]
                # find ORF6 muts
                if "ORF6" in parent.branch_attrs["mutations"]:
                    parents_orf6_muts+=parent.branch_attrs["mutations"]["ORF6"]
                # find ORF7a muts
                if "ORF7a" in parent.branch_attrs["mutations"]:
                    parents_orf7a_muts+=parent.branch_attrs["mutations"]["ORF7a"]
                # find ORF7b muts
                if "ORF7b" in parent.branch_attrs["mutations"]:
                    parents_orf7b_muts+=parent.branch_attrs["mutations"]["ORF7b"]
                # find ORF8 muts
                if "ORF8" in parent.branch_attrs["mutations"]:
                    parents_orf8_muts+=parent.branch_attrs["mutations"]["ORF8"]
                # find ORF9b muts
                if "ORF9b" in parent.branch_attrs["mutations"]:
                    parents_orf9b_muts+=parent.branch_attrs["mutations"]["ORF9b"]
            if hasattr(parent, 'node_attrs') and 'syn_muts' in parent.node_attrs:
                if 'S1' in parent.node_attrs['syn_muts']:
                    parents_s1_syn += parent.node_attrs['syn_muts']['S1']


        # remove reversion mutations from each list
        parents_spike_muts = remove_reversions(parents_spike_muts)
        parents_s1_muts = remove_reversions(parents_s1_muts)
        parents_s2_muts = remove_reversions(parents_s2_muts)
        parents_rdrp_muts = remove_reversions(parents_rdrp_muts)
        parents_n_muts = remove_reversions(parents_n_muts)
        parents_e_muts = remove_reversions(parents_e_muts)
        parents_m_muts = remove_reversions(parents_m_muts)
        parents_orf3a_muts = remove_reversions(parents_orf3a_muts)
        parents_orf6_muts = remove_reversions(parents_orf6_muts)
        parents_orf7a_muts = remove_reversions(parents_orf7a_muts)
        parents_orf7b_muts = remove_reversions(parents_orf7b_muts)
        parents_orf8_muts = remove_reversions(parents_orf8_muts)
        parents_orf9b_muts = remove_reversions(parents_orf9b_muts)

        parents_nsp1_muts = remove_reversions(parents_nsp1_muts)
        parents_nsp2_muts = remove_reversions(parents_nsp2_muts)
        parents_nsp3_muts = remove_reversions(parents_nsp3_muts)
        parents_nsp4_muts = remove_reversions(parents_nsp4_muts)
        parents_nsp5_muts = remove_reversions(parents_nsp5_muts)
        parents_nsp6_muts = remove_reversions(parents_nsp6_muts)
        parents_nsp7_muts = remove_reversions(parents_nsp7_muts)
        parents_nsp8_muts = remove_reversions(parents_nsp8_muts)
        parents_nsp9_muts = remove_reversions(parents_nsp9_muts)
        parents_nsp10_muts = remove_reversions(parents_nsp10_muts)
        parents_nsp13_muts = remove_reversions(parents_nsp13_muts)
        parents_nsp14_muts = remove_reversions(parents_nsp14_muts)
        parents_nsp15_muts = remove_reversions(parents_nsp15_muts)
        parents_nsp16_muts = remove_reversions(parents_nsp16_muts)


        # count deletion of adjacent nucleotides as one mutation event
        spike_mutation_count = consolidate_deletions(parents_spike_muts)
        s1_mutation_count = consolidate_deletions(parents_s1_muts)
        s2_mutation_count = consolidate_deletions(parents_s2_muts)
        rdrp_mutation_count = consolidate_deletions(parents_rdrp_muts)
        n_mutation_count = consolidate_deletions(parents_n_muts)
        e_mutation_count = consolidate_deletions(parents_e_muts)
        m_mutation_count = consolidate_deletions(parents_m_muts)
        orf3a_mutation_count = consolidate_deletions(parents_orf3a_muts)
        orf6_mutation_count = consolidate_deletions(parents_orf6_muts)
        orf7a_mutation_count = consolidate_deletions(parents_orf7a_muts)
        orf7b_mutation_count = consolidate_deletions(parents_orf7b_muts)
        orf8_mutation_count = consolidate_deletions(parents_orf8_muts)
        orf9b_mutation_count = consolidate_deletions(parents_orf9b_muts)

        nsp1_mutation_count = consolidate_deletions(parents_nsp1_muts)
        nsp2_mutation_count = consolidate_deletions(parents_nsp2_muts)
        nsp3_mutation_count = consolidate_deletions(parents_nsp3_muts)
        nsp4_mutation_count = consolidate_deletions(parents_nsp4_muts)
        nsp5_mutation_count = consolidate_deletions(parents_nsp5_muts)
        nsp6_mutation_count = consolidate_deletions(parents_nsp6_muts)
        nsp7_mutation_count = consolidate_deletions(parents_nsp7_muts)
        nsp8_mutation_count = consolidate_deletions(parents_nsp8_muts)
        nsp9_mutation_count = consolidate_deletions(parents_nsp9_muts)
        nsp10_mutation_count = consolidate_deletions(parents_nsp10_muts)
        nsp13_mutation_count = consolidate_deletions(parents_nsp13_muts)
        nsp14_mutation_count = consolidate_deletions(parents_nsp14_muts)
        nsp15_mutation_count = consolidate_deletions(parents_nsp15_muts)
        nsp16_mutation_count = consolidate_deletions(parents_nsp16_muts)


        node.node_attrs["s1_syn_accumulation"] = len(parents_s1_syn)

        node.node_attrs["spike_accumulation"] = spike_mutation_count
        node.node_attrs["s1_accumulation"] = s1_mutation_count
        node.node_attrs["s2_accumulation"] = s2_mutation_count
        node.node_attrs["rdrp_accumulation"] = rdrp_mutation_count
        node.node_attrs["n_accumulation"] = n_mutation_count
        node.node_attrs["e_accumulation"] = e_mutation_count
        node.node_attrs["m_accumulation"] = m_mutation_count
        node.node_attrs["orf3a_accumulation"] = orf3a_mutation_count
        node.node_attrs["orf6_accumulation"] = orf6_mutation_count
        node.node_attrs["orf7a_accumulation"] = orf7a_mutation_count
        node.node_attrs["orf7b_accumulation"] = orf7b_mutation_count
        node.node_attrs["orf8_accumulation"] = orf8_mutation_count
        node.node_attrs["orf9b_accumulation"] = orf9b_mutation_count
        node.node_attrs["nsp1_accumulation"] = nsp1_mutation_count
        node.node_attrs["nsp2_accumulation"] = nsp2_mutation_count
        node.node_attrs["nsp3_accumulation"] = nsp3_mutation_count
        node.node_attrs["nsp4_accumulation"] = nsp4_mutation_count
        node.node_attrs["nsp5_accumulation"] = nsp5_mutation_count
        node.node_attrs["nsp6_accumulation"] = nsp6_mutation_count
        node.node_attrs["nsp7_accumulation"] = nsp7_mutation_count
        node.node_attrs["nsp8_accumulation"] = nsp8_mutation_count
        node.node_attrs["nsp9_accumulation"] = nsp9_mutation_count
        node.node_attrs["nsp10_accumulation"] = nsp10_mutation_count
        node.node_attrs["nsp13_accumulation"] = nsp13_mutation_count
        node.node_attrs["nsp14_accumulation"] = nsp14_mutation_count
        node.node_attrs["nsp15_accumulation"] = nsp15_mutation_count
        node.node_attrs["nsp16_accumulation"] = nsp16_mutation_count
    return tree


def add_del_accumulation_attr(tree):
    """
    For each node on the tree, count the number of deletions that have happened along the path from root to node (including this node). Separate the counts by gene and store them as an attribute of the node. A deletion of several adjacent amino acids is counted as one mutation event.

    Also store an attribute that lists the number of nonsynonymous SNPs accumulated at each node. This is different than the attribute added by `add_mut_accumulation_attr`, which lists a count of SNPs+deletions

    Must run `add_mut_accumulation_attr` on the tree before running `add_del_accumulation_attr`
    """


    for node in tree.find_clades():

        #Find all parents of the node
        parents = get_parent(tree, node)

        #Find mutations that occur in the parents
        parents_spike_dels = []
        parents_s1_dels = []
        parents_s2_dels = []
        parents_rdrp_dels = []
        parents_nsp6_dels = []
        parents_nsp4_dels = []
        parents_n_dels = []
        parents_e_dels = []
        parents_m_dels = []
        parents_orf7a_dels = []

        for parent in parents:
            if hasattr(parent, "branch_attrs") and "mutations" in parent.branch_attrs:
                if "S" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["S"]:
                        if mut[-1] == '-':
                            parents_spike_dels+=[mut]
                            # nextstrain calls pos 13 in sigpep, not S1
                            if int(mut[1:-1]) in range(14,686):
                                parents_s1_dels+=[mut]
                            elif int(mut[1:-1]) in range(687,1274):
                                parents_s2_dels+=[mut]
                #find RdRp muts
                #and Nsp4 and 6 muts
                if "ORF1a" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["ORF1a"]:
                        if mut[-1] == '-':
                            if int(mut[1:-1]) in range(4393,4401):
                                #renumber mut according to rdrp protein
                                rdrp_mut = f'{mut[0]}{int(mut[1:-1])-4393}{mut[-1]}'
                                parents_rdrp_dels+=[rdrp_mut]
                            elif int(mut[1:-1]) in range(3570,3859):
                                # exclude this ancestral mut
                                if mut!= 'K3833N':
                                    #renumber mut according to nsp6 protein
                                    nsp6_mut = f'{mut[0]}{int(mut[1:-1])-3570}{mut[-1]}'
                                    parents_nsp6_dels+=[nsp6_mut]
                            elif int(mut[1:-1]) in range(2754,3263):
                                #renumber mut according to nsp6 protein
                                nsp4_mut = f'{mut[0]}{int(mut[1:-1])-2777}{mut[-1]}'
                                parents_nsp4_dels+=[nsp4_mut]
                if "ORF1b" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["ORF1b"]:
                        if mut[-1] == '-':
                            if int(mut[1:-1]) in range(1,923):
                                #renumber mut according to rdrp protein
                                rdrp_mut = f'{mut[0]}{int(mut[1:-1])+9}{mut[-1]}'
                                parents_rdrp_dels+=[rdrp_mut]
                # find N muts
                if "N" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["N"]:
                        if mut[-1] == '-':
                            parents_n_dels+=parent.branch_attrs["mutations"]["N"]
                # find E muts
                if "E" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["E"]:
                        if mut[-1] == '-':
                            parents_e_dels+=parent.branch_attrs["mutations"]["E"]
                # find M muts
                if "M" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["M"]:
                        if mut[-1] == '-':
                            parents_m_dels+=parent.branch_attrs["mutations"]["M"]
                # find ORF7a muts
                if "ORF7a" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["ORF7a"]:
                        if mut[-1] == '-':
                            parents_m_dels+=parent.branch_attrs["mutations"]["ORF7a"]



        # count deletion of adjacent nucleotides as one mutation event
        spike_deletion_count = consolidate_deletions(parents_spike_dels)
        s1_deletion_count = consolidate_deletions(parents_s1_dels)
        s2_deletion_count = consolidate_deletions(parents_s2_dels)
        rdrp_deletion_count = consolidate_deletions(parents_rdrp_dels)
        nsp6_deletion_count = consolidate_deletions(parents_nsp6_dels)
        nsp4_deletion_count = consolidate_deletions(parents_nsp4_dels)
        n_deletion_count = consolidate_deletions(parents_n_dels)
        e_deletion_count = consolidate_deletions(parents_e_dels)
        m_deletion_count = consolidate_deletions(parents_m_dels)
        orf7a_deletion_count = consolidate_deletions(parents_orf7a_dels)

        # add node attribute that lists the number of deletions accumulated at each node
        node.node_attrs['deletion_accumulation']= {'Spike': spike_deletion_count,
                                                   'S1': s1_deletion_count, 'S2':s2_deletion_count,
                                                   'RdRp':rdrp_deletion_count, 'Nsp6':nsp6_deletion_count,
                                                   'Nsp4':nsp4_deletion_count, 'N':n_deletion_count,
                                                   'M': m_deletion_count, 'E': e_deletion_count,
                                                   'ORF7a': orf7a_deletion_count}
        # add node attribute that lists the number of nonsynonymous SNPs accumulated at each node
        node.node_attrs['nonsyn_snps_accumulation'] = {'Spike': node.node_attrs["spike_accumulation"] - spike_deletion_count,
                                                       'S1': node.node_attrs["s1_accumulation"] - s1_deletion_count,
                                                       'S2': node.node_attrs["s2_accumulation"] - s2_deletion_count,
                                                       'RdRp': node.node_attrs["rdrp_accumulation"] - rdrp_deletion_count,
                                                       'Nsp6': node.node_attrs["nsp6_accumulation"] - nsp6_deletion_count,
                                                       'Nsp4': node.node_attrs["nsp4_accumulation"] - nsp4_deletion_count,
                                                       'N': node.node_attrs["n_accumulation"] - n_deletion_count,
                                                       'M': node.node_attrs["m_accumulation"] - m_deletion_count,
                                                       'E': node.node_attrs["e_accumulation"] - e_deletion_count,
                                                       'ORF7a': node.node_attrs["orf7a_accumulation"] - orf7a_deletion_count}

    return tree

def add_mut_at_node_attr(tree):
    """
    For each node, find the number of mutations that happened within each gene. Store this as an attribute of the node
    Nsp genome locations from: https://www.ncbi.nlm.nih.gov/gene/43740578
    """
    for node in tree.find_clades(terminal=False):

        node.nonsyn_at_node = {}
        node.syn_at_node = {}

        s1_syn_at_this_node = []
        if hasattr(node, "node_attrs") and 'S1' in node.node_attrs['syn_muts']:
            s1_syn_at_this_node.append(node.node_attrs['syn_muts']['S1'])
        node.syn_at_node['S1'] = len(s1_syn_at_this_node)


        if hasattr(node, 'branch_attrs'):

            s1_nonsyn_at_this_node = []
            s2_nonsyn_at_this_node = []
            if "S" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["S"]:
                    if int(mut[1:-1]) in range(14,686):
                        s1_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(687,1274):
                        s2_nonsyn_at_this_node.append(mut)

            s1_consolidated = consolidate_deletions_2(s1_nonsyn_at_this_node)
            node.nonsyn_at_node['S1'] = len(s1_consolidated)
            s2_consolidated = consolidate_deletions_2(s2_nonsyn_at_this_node)
            node.nonsyn_at_node['S2'] = len(s2_consolidated)



            rdrp_nonsyn_at_this_node = []

            nsp1_nonsyn_at_this_node = []
            nsp2_nonsyn_at_this_node = []
            nsp3_nonsyn_at_this_node = []
            nsp4_nonsyn_at_this_node = []
            nsp6_nonsyn_at_this_node = []
            nsp5_nonsyn_at_this_node = []
            nsp7_nonsyn_at_this_node = []
            nsp8_nonsyn_at_this_node = []
            nsp9_nonsyn_at_this_node = []
            nsp10_nonsyn_at_this_node = []
            nsp13_nonsyn_at_this_node = []
            nsp14_nonsyn_at_this_node = []
            nsp15_nonsyn_at_this_node = []
            nsp16_nonsyn_at_this_node = []
            if "ORF1a" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["ORF1a"]:
                    if int(mut[1:-1]) in range(4393,4401):
                        rdrp_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(1,180):
                        nsp1_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(181,818):
                        nsp2_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(819,2763):
                        nsp3_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(2764,3263):
                        nsp4_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(3264,3569):
                        nsp5_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(3570,3859):
                        # exclude this ancestral mut
                        if mut!= 'K3833N':
                            nsp6_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(3860,3942):
                        nsp7_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(3943,4140):
                        nsp8_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(4141,4253):
                        nsp9_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(4254,4392):
                        nsp10_nonsyn_at_this_node.append(mut)

            if "ORF1b" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["ORF1b"]:
                    if int(mut[1:-1]) in range(1,923):
                        rdrp_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(924,1524):
                        nsp13_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(1525,2051):
                        nsp14_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(2052,2397):
                        nsp15_nonsyn_at_this_node.append(mut)
                    elif int(mut[1:-1]) in range(2398,2695):
                        nsp16_nonsyn_at_this_node.append(mut)

            rdrp_consolidated = consolidate_deletions_2(rdrp_nonsyn_at_this_node)
            node.nonsyn_at_node['RdRp'] =  len(rdrp_consolidated)

            nsp1_consolidated = consolidate_deletions_2(nsp1_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp1'] = len(nsp1_consolidated)
            nsp2_consolidated = consolidate_deletions_2(nsp2_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp2'] = len(nsp2_consolidated)
            nsp3_consolidated = consolidate_deletions_2(nsp3_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp3'] = len(nsp3_consolidated)
            nsp4_consolidated = consolidate_deletions_2(nsp4_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp4'] = len(nsp4_consolidated)
            nsp5_consolidated = consolidate_deletions_2(nsp5_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp5'] = len(nsp5_consolidated)
            nsp6_consolidated = consolidate_deletions_2(nsp6_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp6'] = len(nsp6_consolidated)
            nsp7_consolidated = consolidate_deletions_2(nsp7_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp7'] = len(nsp7_consolidated)
            nsp8_consolidated = consolidate_deletions_2(nsp8_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp8'] = len(nsp8_consolidated)
            nsp9_consolidated = consolidate_deletions_2(nsp9_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp9'] = len(nsp9_consolidated)
            nsp10_consolidated = consolidate_deletions_2(nsp10_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp10'] = len(nsp10_consolidated)
            nsp13_consolidated = consolidate_deletions_2(nsp13_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp13'] = len(nsp13_consolidated)
            nsp14_consolidated = consolidate_deletions_2(nsp14_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp14'] = len(nsp14_consolidated)
            nsp15_consolidated = consolidate_deletions_2(nsp15_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp15'] = len(nsp15_consolidated)
            nsp16_consolidated = consolidate_deletions_2(nsp16_nonsyn_at_this_node)
            node.nonsyn_at_node['Nsp16'] = len(nsp16_consolidated)

            n_nonsyn_at_this_node = []
            if "N" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["N"]:
                    n_nonsyn_at_this_node.append(mut)

            n_consolidated = consolidate_deletions_2(n_nonsyn_at_this_node)
            node.nonsyn_at_node['N'] = len(n_consolidated)

            e_nonsyn_at_this_node = []
            if "E" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["E"]:
                    e_nonsyn_at_this_node.append(mut)

            e_consolidated = consolidate_deletions_2(e_nonsyn_at_this_node)
            node.nonsyn_at_node['E'] = len(e_consolidated)

            m_nonsyn_at_this_node = []
            if "M" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["M"]:
                    m_nonsyn_at_this_node.append(mut)

            m_consolidated = consolidate_deletions_2(m_nonsyn_at_this_node)
            node.nonsyn_at_node['M'] = len(m_consolidated)

            orf3a_nonsyn_at_this_node = []
            if "ORF3a" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["ORF3a"]:
                    orf3a_nonsyn_at_this_node.append(mut)

            orf3a_consolidated = consolidate_deletions_2(orf3a_nonsyn_at_this_node)
            node.nonsyn_at_node['ORF3a'] = len(orf3a_consolidated)

            orf6_nonsyn_at_this_node = []
            if "ORF6" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["ORF6"]:
                    orf6_nonsyn_at_this_node.append(mut)

            orf6_consolidated = consolidate_deletions_2(orf6_nonsyn_at_this_node)
            node.nonsyn_at_node['ORF6'] = len(orf6_consolidated)

            orf7a_nonsyn_at_this_node = []
            if "ORF7a" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["ORF7a"]:
                    orf7a_nonsyn_at_this_node.append(mut)

            orf7a_consolidated = consolidate_deletions_2(orf7a_nonsyn_at_this_node)
            node.nonsyn_at_node['ORF7a'] = len(orf7a_consolidated)

            orf7b_nonsyn_at_this_node = []
            if "ORF7b" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["ORF7b"]:
                    orf7b_nonsyn_at_this_node.append(mut)

            orf7b_consolidated = consolidate_deletions_2(orf7b_nonsyn_at_this_node)
            node.nonsyn_at_node['ORF7b'] = len(orf7b_consolidated)

            orf8_nonsyn_at_this_node = []
            if "ORF8" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["ORF8"]:
                    orf8_nonsyn_at_this_node.append(mut)

            orf8_consolidated = consolidate_deletions_2(orf8_nonsyn_at_this_node)
            node.nonsyn_at_node['ORF8'] = len(orf8_consolidated)

            orf9b_nonsyn_at_this_node = []
            if "ORF9b" in node.branch_attrs["mutations"]:
                for mut in node.branch_attrs["mutations"]["ORF9b"]:
                    orf9b_nonsyn_at_this_node.append(mut)

            orf9b_consolidated = consolidate_deletions_2(orf9b_nonsyn_at_this_node)
            node.nonsyn_at_node['ORF9b'] = len(orf9b_consolidated)
    return tree


def add_changes_from_ref_attr(tree):
    """
    For each node, find all mutations that occurred between the root and the node. Separate these mutations by gene and by synonymous/nonsynonymous. Add a dictionary of these mutations as an attribute to each node
    """
    for node in tree.find_clades():
        #Find all parents of the node (includes node too)
        parents = get_parent(tree, node)

        #Find mutations that occur in the parents
        parents_spike_muts = []
        parents_s1_muts = []
        parents_s2_muts = []
        parents_rdrp_muts = []
        parents_n_muts = []
        parents_e_muts = []
        parents_m_muts = []

        parents_s1_syn = []
        parents_s2_syn = []
        parents_rdrp_syn = []
        parents_spike_syn = []
        parents_n_syn = []
        parents_e_syn = []
        parents_m_syn = []


        for parent in parents:
            if hasattr(parent, "branch_attrs") and "mutations" in parent.branch_attrs:
                if "S" in parent.branch_attrs["mutations"]:
                    parents_spike_muts+=parent.branch_attrs["mutations"]["S"]
                    for mut in parent.branch_attrs["mutations"]["S"]:
                        # nextstrain calls pos 13 in sigpep, not S1
                        if int(mut[1:-1]) in range(14,686):
                            parents_s1_muts+=[mut]
                        elif int(mut[1:-1]) in range(687,1274):
                            parents_s2_muts+=[mut]
                #find RdRp muts
                #and Nsp4 and 6 muts
                if "ORF1a" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["ORF1a"]:
                        if int(mut[1:-1]) in range(4393,4401):
                            #renumber mut according to rdrp protein
                            rdrp_mut = f'{mut[0]}{int(mut[1:-1])-4492}{mut[-1]}'
                            parents_rdrp_muts+=[rdrp_mut]
                if "ORF1b" in parent.branch_attrs["mutations"]:
                    for mut in parent.branch_attrs["mutations"]["ORF1b"]:
                        if int(mut[1:-1]) in range(1,923):
                            #renumber mut according to rdrp protein
                            rdrp_mut = f'{mut[0]}{int(mut[1:-1])+9}{mut[-1]}'
                            parents_rdrp_muts+=[rdrp_mut]
                # find N muts
                if "N" in parent.branch_attrs["mutations"]:
                    parents_n_muts+=parent.branch_attrs["mutations"]["N"]
                # find E muts
                if "E" in parent.branch_attrs["mutations"]:
                    parents_e_muts+=parent.branch_attrs["mutations"]["E"]
                # find M muts
                if "M" in parent.branch_attrs["mutations"]:
                    parents_m_muts+=parent.branch_attrs["mutations"]["M"]
            if hasattr(parent, 'node_attrs') and 'syn_muts' in parent.node_attrs:
                if 'S' in parent.node_attrs['syn_muts']:
                    parents_spike_syn += parent.node_attrs['syn_muts']['S']
                if 'S1' in parent.node_attrs['syn_muts']:
                    parents_spike_syn += parent.node_attrs['syn_muts']['S1']
                    parents_s1_syn += parent.node_attrs['syn_muts']['S1']
                if 'S2' in parent.node_attrs['syn_muts']:
                    parents_spike_syn += parent.node_attrs['syn_muts']['S2']
                    parents_s2_syn += parent.node_attrs['syn_muts']['S2']
                if 'RdRp' in parent.node_attrs['syn_muts']:
                    parents_rdrp_syn += parent.node_attrs['syn_muts']['RdRp']
                if 'E' in parent.node_attrs['syn_muts']:
                    parents_e_syn += parent.node_attrs['syn_muts']['E']
                if 'N' in parent.node_attrs['syn_muts']:
                    parents_n_syn += parent.node_attrs['syn_muts']['N']
                if 'M' in parent.node_attrs['syn_muts']:
                    parents_m_syn += parent.node_attrs['syn_muts']['M']


        # remove reversion mutations from each list
        parents_spike_muts = remove_reversions(parents_spike_muts)
        parents_s1_muts = remove_reversions(parents_s1_muts)
        parents_s2_muts = remove_reversions(parents_s2_muts)
        parents_rdrp_muts = remove_reversions(parents_rdrp_muts)
        parents_n_muts = remove_reversions(parents_n_muts)
        parents_e_muts = remove_reversions(parents_e_muts)
        parents_m_muts = remove_reversions(parents_m_muts)
        parents_s1_syn = remove_reversions(parents_s1_syn)
        parents_s2_syn = remove_reversions(parents_s2_syn)
        parents_rdrp_syn = remove_reversions(parents_rdrp_syn)
        parents_spike_syn = remove_reversions(parents_spike_syn)
        parents_e_syn = remove_reversions(parents_e_syn)
        parents_n_syn = remove_reversions(parents_n_syn)
        parents_m_syn = remove_reversions(parents_m_syn)


        # count deletion of adjacent nucleotides as one mutation event
        spike_mutation_list = consolidate_deletions_2(parents_spike_muts)
        s1_mutation_list = consolidate_deletions_2(parents_s1_muts)
        s2_mutation_list = consolidate_deletions_2(parents_s2_muts)
        rdrp_mutation_list = consolidate_deletions_2(parents_rdrp_muts)
        n_mutation_list = consolidate_deletions_2(parents_n_muts)
        e_mutation_list = consolidate_deletions_2(parents_e_muts)
        m_mutation_list = consolidate_deletions_2(parents_m_muts)
        s1_syn_mutation_list = consolidate_deletions_2(parents_s1_syn)
        s2_syn_mutation_list = consolidate_deletions_2(parents_s2_syn)
        rdrp_syn_mutation_list = consolidate_deletions_2(parents_rdrp_syn)
        spike_syn_mutation_list = consolidate_deletions_2(parents_spike_syn)
        n_syn_mutation_list = consolidate_deletions_2(parents_n_syn)
        e_syn_mutation_list = consolidate_deletions_2(parents_e_syn)
        m_syn_mutation_list = consolidate_deletions_2(parents_m_syn)

        node.node_attrs["changes_from_ref"] = {'s1_non':s1_mutation_list, 's2_non':s2_mutation_list, 'rdrp_non':rdrp_mutation_list, 'spike_non':spike_mutation_list, 'e_non':e_mutation_list, 'n_non':n_mutation_list, 'm_non':m_mutation_list, 's1_syn': s1_syn_mutation_list, 's2_syn':s2_syn_mutation_list, 'rdrp_syn': rdrp_syn_mutation_list, 'spike_syn': spike_syn_mutation_list, 'e_syn': e_syn_mutation_list, 'n_syn':n_syn_mutation_list, 'm_syn': m_syn_mutation_list}
    return tree


def final_pos_genotype(mutation_list, mutation):
    """
    Given a list of mutations that occurred at certain position,
    find the whether the final genotype at this position is the specified mutation or not (return True or False)
    """

    if len(mutation_list) == 0:
        mutation_at_pos = False

    elif len(mutation_list) == 1:
        if mutation_list[0][-1] == mutation:
            mutation_at_pos = True
        else:
            mutation_at_pos = False


    elif len(mutation_list) > 1:
        if mutation_list[-1][-1] == mutation:
            mutation_at_pos = True
        else:
            mutation_at_pos = False

    return mutation_at_pos


def prune_tree(tree, num_descending_tips):
    """
    Return a list of all terminal branches on a tree that is
    pruned to only branches with 3 or more descending tips
    """

    # make a list of all terminal nodes on a tree truncated to branches
    # with at least num_descending_tips descending tips
    pruned_terminals = []
    for node in tree.find_clades(terminal=False):
        # make sure node has at least the specified number of descending tips
        if len(node.get_terminals()) >= num_descending_tips:

            # if a parent of the node is already listed, overwrite it
            parents = get_parent(tree, node)[:-1]
            parent_names = [p.name for p in parents]

            # want path ends, so parents true terminals should not be in pruned_terminals
            # find whether any parents are listed in pruned_terminals
            matches = list((i in pruned_terminals for i in parent_names))

            # if a parent of this node is in pruned_terminals,
            # that parent is not a true terminal and should be removed
            if True in matches:
                match_indicies = [index for index, value in enumerate(matches) if value == True]
                for match_index in match_indicies:
                    pruned_terminals.remove(parent_names[match_index])

            pruned_terminals.append(node.name)

    # remove root
    for x in ['NODE_0000000', 'NODE_0000001']:
        if x in pruned_terminals:
            pruned_terminals.remove(x)


    return pruned_terminals

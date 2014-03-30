# -*- coding: utf-8 -*-
from collections import defaultdict, Counter
import itertools
import logging
import sys

__author__ = "Sergey Aganezov"
__email__ = "aganezov@gwu.edu"
__status__ = "develop"

"""
This script shall be invoked as a command line tool.
This script is written with python3.4+ support.
This script accepts only command line arguments.

Idea:
    the main idea goes as follows. As gene family mapping has to be consistent among different phylogenetic tree
    splits, we can double check that, taking several genome from more recent split in the tree, and more general
    genome set from more dated split in the tree. Take two orthology mapping files, and check same gene families
    along with particular genome genes (coding exons) that are mapped to them. Inconsistency arises in the case,
    when we have same gene (coding exons id) from particular genome is mapped to a different gene families in
    different orthology map.

    As in different mapping files gene families names might differ, further remapping gene families must be performed.
    The following scheme for this is executed:
        1. for the first mapping file map all gene families, sorted in alphabetic order, to respective integers
        2. for the second file, get organism, that is present in both original and observed file, for that organism,
        using the gene_id column, map gene families into each other, and then into respective integers.
        3. after that for the rest of organism, using created in step 2 mapping, check similarity if gene_id mapping
        into gene families
        3.1 if for some organism, there if gene family in file 2, that wasn't mapped at step 2, try to find such
        organism in file one, try to find respective gene_id, and perform mapping, similar to step 2


Input:
    a list of source files with orthology mapping supplied with space delimiter as cmd arguments for this script
    list shall contain orthology mapping data in tabtext format
    example:
        >>> python3.4 orthology_mapping_file_1 orthology_mapping_file_2 ...

Output:
    the script writes to the standard output
    the script outputs only the inconsistent data
"""


# is used for assigning new integers mapping values to gene families
integers_source = iter(itertools.count())


def retrieve_info_from_string(raw_data_string):
    """
    retrieves relative information from raw tab separated string from orthology mapping file
    Args:
        raw_data_string: original tab-separated string

    Returns:
        tuple of following columns from original string: (gene family id, gene id, organism name)

    Raises:
        IndexError, as splitting might not go flawless
    """
    data = raw_data_string.split("\t")
    return data[1], data[3], data[4]


def split_file_by_o_ids(file_name):
    """ Parses original file, retrieves relative information, that splits overall data according to orthology ids
    Args:
        file_name:

    Returns:
            defaultdict: {gene family id: [list of gene ids, organisms for this particular gene family id]}
    """
    with open(file_name, "r") as source:
        o_id_splitting = defaultdict(list)
        data = source.readlines()
        retrieved_data = []
        for cnt, line in enumerate(data):
            if cnt > 0:
                try:
                    retrieved_data.append(retrieve_info_from_string(line))
                except IndexError:
                    continue
        for o_id, g_id, organism in retrieved_data:
            o_id_splitting[o_id].append((g_id, organism))
        return o_id_splitting


def main(source_files_list, dest=None):
    gene_family_id_mapping = {}
    gene_id_mapping = defaultdict(dict)

    retrieve_o_id_value = lambda o_id, storage: storage.get(o_id, None)
    retrieve_g_id_value = lambda g_id, organism, storage: storage[organism].get(g_id, None)
    retrieve_first_value = lambda iterable: next(filter(lambda x: x is not None, iterable))

    for source_file in source_files_list:
        o_id_split_data = split_file_by_o_ids(source_file)

        miss_matched_gene_families_cnt = 0
        miss_matched_gene_families_freq = defaultdict(list)
        number_of_genes_per_family = {}

        for o_id in sorted(o_id_split_data.keys()):  # sorting is done for data output consistency
            g_id_values = list(
                retrieve_g_id_value(g_id, organism, gene_id_mapping) for g_id, organism in o_id_split_data[o_id])
            number_of_genes_per_family[o_id] = len(g_id_values)
            g_id_values_cnt = Counter(g_id_value for g_id_value in g_id_values if g_id_value is not None)
            o_id_value = retrieve_o_id_value(o_id, gene_family_id_mapping)
            #########################################################################################################
            # in this case we have the situation when neither gene family was previously mapped
            # and non of the gene ids, that are mapped to the current gene family, were previously mapped
            #########################################################################################################
            if o_id_value is None and all(g_id is None for g_id in g_id_values):

                logging.debug("o_id {o_id} and none of g_ids {g_ids} were previously mapped".format(
                    o_id=o_id,
                    g_ids="(" + ", ".join(str(g_id) for g_id, organism in o_id_split_data[o_id]) + ")"))

                next_int = next(integers_source)
                gene_family_id_mapping[o_id] = next_int
                for g_id, organism in o_id_split_data[o_id]:
                    gene_id_mapping[organism][g_id] = next_int
                continue
            #########################################################################################################
            # in this case we have a situation, when gene family was not previously mapped, but
            # some of gene ids, that are mapped to the current gene family, were previously mapped
            # used first gene id mapping as mapping for gene family
            #########################################################################################################
            if o_id_value is None and not all(g_id_value is None for g_id_value in g_id_values):

                logging.debug("o_id {o_id} was not previously mapped, but some of g_ids were".format(o_id=o_id))

                o_id_value = g_id_values_cnt.most_common()[0][0]  # extract the actual value, of the most common key
                gene_family_id_mapping[o_id] = o_id_value
            #########################################################################################################
            # at this point gene family is mapped to some value for sure, and we need to check, if gene id mapping
            # is consistent with gene family mapping
            #########################################################################################################
            for g_id_value, (g_id, organism) in zip(g_id_values, o_id_split_data[o_id]):
                if g_id_value is None:
                    gene_id_mapping[organism][g_id] = o_id_value
                elif o_id_value != g_id_value:
                    pass
                    # print("Miss matching!!! gene_id {g_id} is miss matched!!".format(g_id=g_id), file=dest)
            #########################################################################################################
            #
            #########################################################################################################
            if len(g_id_values_cnt) > 1:
                key = sum(x[1] for x in g_id_values_cnt.most_common()[1:])  # number of genes, for g_ids, different
                                                                            # from the most popular
                miss_matched_gene_families_cnt += 1
                miss_matched_gene_families_freq[key].append(o_id)
        print("Orthology mapping file {file_name} contained {mgfc} gene families (out of {ogfc}), where at least"
              " one gene was mapped, differently, from previously observed orthology mapping files."
              "".format(file_name=source_file, mgfc=miss_matched_gene_families_cnt,
                        ogfc=len(number_of_genes_per_family)),
              file=dest)
        if len(miss_matched_gene_families_freq) > 0:
            print("Among miss-mapped gene families, there were", file=dest)
            for key, value in sorted(miss_matched_gene_families_freq.items(), key=lambda item: -item[0]):
                print("\t{value} gene families have {key} miss-mapped genes. List of these families"
                      "".format(value=len(value), key=key), file=dest)
                tmp = map(lambda x: "{o_id} (out of {overall})".format(o_id=x, overall=number_of_genes_per_family[x]),
                          value)
                print("\t\t" + "\n\t\t".join(str_value for str_value in tmp), file=dest)



if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    orthology_source_files = cmd_args
    logging.basicConfig(level=50)
    main(orthology_source_files)




#!/usr/local/bin/python3.4
# -*- coding: utf-8 -*-
from collections import defaultdict, Counter
import itertools
import sys

__author__ = "Sergey Aganezov"
__email__ = "aganezov@gwu.edu"
__status__ = "develop"

""" Overview / help

This script shall be invoked as a command line tool.
This script is written with python3.4+ support.
This script accepts only command line arguments.
This script can be used in a Unix/Linux pipe, as it produces simple output, but must be invoked as a first producer.

Idea:
    the main idea goes as follows. As gene family mapping has to be consistent among different phylogenetic tree
    splits, we can double check that, taking several genomes from more recent split in the tree, and more general
    genome set from more dated split in the tree. Take two orthology mapping files, and check same gene families
    along with particular genome genes (coding exons) that are mapped to them. Inconsistency arises in the case,
    when we have same gene (coding exons id) from particular genome mapped to a different gene families in
    different orthology map.

    As in different mapping files gene families names might differ, further remapping gene families must be performed.
    The following scheme for this is executed:
        1. for the first mapping file we map all gene families, sorted in alphabetic order, to respective integers
        2. for the second file we split data by gene families id, and then we try to determine, to which integer
        this gene family must be mapped. This value is determined by the gene, that belong to this family (if any
        of them were previously mapped). After that we double check for all gene, that belong to this particular family,
        if they were mapped previously, and if so, if they were mapped to the same integer, this gene family, as a whole
        if mapped.

    Example:
        file 1 >>>
            gene_family_id_1 gene_1
            gene_family_id_1 gene_2
            gene_family_id_1 gene_3
            gene_family_id_1 gene_4
            gene_family_id_1 gene_5
            gene_family_id_2 gene_11
            gene_family_id_2 gene_21
            gene_family_id_2 gene_31
            gene_family_id_2 gene_41
            gene_family_id_2 gene_51

        file >>>
            gene_family_id_11 gene_1
            gene_family_id_11 gene_2
            gene_family_id_11 gene_3
            gene_family_id_11 gene_4
            gene_family_id_11 gene_51
            gene_family_id_21 gene_5
            gene_family_id_21 gene_11
            gene_family_id_21 gene_21
            gene_family_id_21 gene_31
            gene_family_id_21 gene_41


        first file gets processed, and gene_family_id_1 gets mapped to 1, as well, as all genes, that belong to it.
        gene_family_id_2 gets mapped to 2, as well as all genes, that belong to it.

        during the processing of the second file, gene_family_id_11 gets mapped to 1, as 4, out of 5 genes, that
        belong to it, were previously mapped to 1. after that each gene in gene_family_id_11 gets checked, to see, if
        it was previously checked, and if so, checked if it was mapped to 1. Fro this process gene_51 gets discovered,
        as it was previously mapped to 2, but here belongs to family, that is mapped to 1. Same analysis goes for
        gene_family_id_21, and gene_5 gets discovered.

Input:
    2 source files with orthology mapping supplied with space delimiter as cmd arguments for this script
    list shall contain orthology mapping data in tabtext format
    example:
        >>> python3.4 check.py orthology_mapping_file_1 orthology_mapping_file_2

Output:
    the script writes to the standard output
    the script outputs only the orthology ids of inconsistent gene families

    Example:
        Orthology mapping file orthology_mapping_file_1 contained 0 gene families (out of 19000), where at least one
        gene was mapped, differently, from previously observed orthology mapping files.

        Orthology mapping file orthology_mapping_file_2 contained 200 gene families (out of 16754), where at least
        one gene was mapped, differently, from previously observed orthology mapping files.
        Among miss-mapped gene families, there were
            2 gene families have 1 miss-mapped genes. List of these families
                gene_family_id_11 (out of 5)
                gene_family_id_21 (out of 5)
"""


# is used for assigning new integers mapping values to gene families
integers_source = iter(itertools.count())


def retrieve_info_from_string(raw_data_string):
    """  Retrieves relative information from raw, tab separated string from orthology mapping file

    At this moment we're only interested in "gene family id", "gene id" and "organism name" columns

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

    starting with row 2, as first row must be dedicated to column naming
    skips line, that were not be properly processed with retrieving information procedure

    Args:
        file_name: source orthology tabtext (tab separated text) file.

    Returns:
            defaultdict: {gene family id: list of tuples (gene ids, organisms name)}
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
            # and non of the gene ids, that belong to the current gene family, were previously mapped
            #########################################################################################################
            if o_id_value is None and all(g_id is None for g_id in g_id_values):
                next_int = next(integers_source)
                gene_family_id_mapping[o_id] = next_int
                for g_id, organism in o_id_split_data[o_id]:
                    gene_id_mapping[organism][g_id] = next_int
                continue
            #########################################################################################################
            # in this case we have a situation, when gene family was not previously mapped, but
            # some of gene ids, that belong to the current gene family, were previously mapped
            # used determine the biggest group of mapped gene ids, and use their value for mapping current gene family
            #########################################################################################################
            if o_id_value is None and not all(g_id_value is None for g_id_value in g_id_values):
                o_id_value = g_id_values_cnt.most_common()[0][0]  # extract the actual value, of the most common key
                gene_family_id_mapping[o_id] = o_id_value
            #########################################################################################################
            # at this point gene family is mapped to some value for sure, we need to mapped all gene ids, that were not
            # previously mapped to the value of the current gene family
            #########################################################################################################
            for g_id_value, (g_id, organism) in zip(g_id_values, o_id_split_data[o_id]):
                if g_id_value is None:
                    gene_id_mapping[organism][g_id] = o_id_value
            #########################################################################################################
            # need to see how many gene ids were miss mapped.
            # as we counted mapping for each gene family, and used the most common one for the gene family mapping
            # all other mapping become inconsistent. count them, to determine how many gene ids are miss mapped in
            # current gene family. And count the number of gene families, which contain at least one miss mapping
            #########################################################################################################
            if len(g_id_values_cnt) > 1:
                key = sum(x[1] for x in g_id_values_cnt.most_common()[1:])  # number of gene ids that are mapped,
                                                                            # differently from the gene family mapping
                miss_matched_gene_families_cnt += 1
                miss_matched_gene_families_freq[key].append(o_id)


        # print("Orthology mapping file {file_name} contained {mgfc} gene families (out of {ogfc}), where at least"
        #       " one gene was mapped, differently, from previously observed orthology mapping files."
        #       "".format(file_name=source_file, mgfc=miss_matched_gene_families_cnt,
        #                 ogfc=len(number_of_genes_per_family)),
        #       file=dest)

        # for each file, that contains at least one gene family, that was miss-mapped, print detailed statistics
        if len(miss_matched_gene_families_freq) > 0:
            # print("Among miss-mapped gene families, there were", file=dest)
            for key, value in sorted(miss_matched_gene_families_freq.items(), key=lambda item: -item[0]):
                # print("\t{value} gene families have {key} miss-mapped genes. List of these families"
                #       "".format(value=len(value), key=key), file=dest)
                tmp = map(lambda x: "{o_id}".format(o_id=x),
                          value)
                print("\n".join(str_value for str_value in tmp), file=dest)


if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    if len(cmd_args) != 2:
        sys.exit(-1)      #  for now program works with only two orthology mapping files
    orthology_source_files = cmd_args
    main(orthology_source_files)




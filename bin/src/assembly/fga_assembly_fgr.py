# -*- coding: utf-8 -*-
from collections import defaultdict
import os
import sys
import re

__author__ = 'Sergey Aganezov'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'

""" Overview / help

This script shall be invoked as a command line tool.
This script is written with python3.4+ support.
This script accepts only command line arguments.
This script cannot be used in a Unix/Linux pipe.

Idea:
    The idea is to take MGRA-log-kind-of file (MGRAL), several gff files, and gene mapping file,
    to report report a genome as a new set of fragment, some of which might have been glued together.

    Script analyses MGRAL file, identifying fusion operations and vertices they affect, after that script maps MGRA
    BG vertices representation with actual gene ids they correspond to. After that gff files are processed in order
    to map respective genes to fragments they are located on, perform several checks, and report which fragments were
    glued together and in which directions. Then each genome gets return as a new set of fragments, some of which are
    composed by a set of glued smaller fragments.

    Example:
        MGRAL file states:
                R=AARA=Arabiensis
                ...
                (1t, oo)x(oo, 2t):{R}

        Gene mapping file states:
                AARA0078 1
                AARA0098 2
        gff file for Arabiensis states:
            K2345 ... AARA0078  .. 1 5 +
            K2356 ... AARA0098  .. 1 10 +

        script will report:
            arabiensis:
                K2345(-) <==> K2356(+)

    double checks:
        script checks that genes, gluing is performed on, are left/right most genes on fragments
        script checks that genes extremities, gluing is performed on, represent the end of the segment
        script check that genes, gluing is performed on, are on different fragments


Input:
    1 MGRAL file
    1 gene mapping file
    n gff files

    MGRAL file format.
    Such file must contain two types of lines:
            1. Y=ABCS=abcdefg
                where Y -- one letter\digit alias for organism name
                where ABCS -- gene_naming prefix, that is used in gff files
                where abcdefg -- genome verbose name
            2. (1h,oo)x(2t,oo):{Y}  (3h,oo)x(4t,oo):{Y}
                where 1h, 2t, oo -- vertices in breakpoint graph
                where Y -- genome one letter/digit alias
                the (a,b)x(c,d):{E} -- gluing operation

    Gene mapping file format
    Such file sut contain only one type of line:
        ABCDEFGD 123
            where ABCDEFGE -- gene id (ABCD here stand for gene_naming prefix, that must match one in MGRAL file)
            where 123 -- value this gene was mapped to, during the transformation to breakpoint graph

    gff files format can be read on website http://www.ensembl.org/info/website/upload/gff.html

Output:

"""

# yeap, evil
# these regular expressions are used to parse MGRAL format file
vrregex = re.compile("\((?:(?P<vertex1>\d+(?:h|t))|oo),(?:(?P<vertex2>\d+(?:h|t))|oo)\)x\((?:(?P<vertex3>\d+(?:h|t))|oo),(?:(?P<vertex4>\d+(?:h|t))|oo)\):\{(?P<genome_alias>\w)\}")
gnsregex = re.compile("(?P<ola>\w)=(?P<gkn>\w+)=(?P<gn>\w+)")


def retrieve_gene_number_mapping_from_file(gene_mapping_file):
    """ Retrieves and separates and returns gene mapping information from supplied file

    Args:
        gene_mapping_file: text file to be parsed. gene_mapping file is specified in Overview docstring

    Returns:
        dict of dicts, where first dict maps 4 letter/digit gene_naming prefix to all genes, that have it,
        and for every such mapping further mapping is performed, where each integer, gene is mapped to, is mapped to
        respective gene (reverse mapping of sorts)

    Raises:
        IndexError: as split might not be performed flawless
    """
    result = defaultdict(dict)
    with open(gene_mapping_file, "r") as source:
        for line in source:
            line = line.strip()
            gene_id, number = line.split()
            number = int(number)
            key = gene_id[:4]
            result[key][number] = gene_id
    return result


def remove_tandem_duplications(fragment_content):
    """ Folds repetitive entries from supplied list to single ones

    at the end we perform a double check, making sure, that in the original sequence we didn't have non-repetitive
    as it shall not be so in this analysis. Though this is just "assert" check which can be disabled.

    Args:
        fragment_content -- list of tuples, where first element represents the identifying value for each entry

    Returns:
        folded version of the list, where each repetitive sequence is folded into one entry

    Raises:
        TypeError: as list might contain not tuples, but some unsubscriptable objects
        AssertionError: if segment contained non-consecutive sequences of genes
    """
    if len(fragment_content) <= 1:
        return fragment_content
    pr_gene = fragment_content[0]
    result = [pr_gene]
    for gene in fragment_content[1:]:
        gene_id = gene[0]
        if gene_id != pr_gene[0]:
            result.append(gene)
            pr_gene = gene
    assert len(set(x[0] for x in result)) == len(result)
    return result


def retrieve_glued_vertices_from_file(mgral_file):
    """ retrieves and processes information from mgral like file

    gathers information about genome naming: for ech genome information about it one letter/digit alias is retrieved,
    and information about its gene_prefix naming is retrieved

    for each genome, statistics regarding gluing operations, that were performed in it is retrieved

    double check is performed to make sure, that glued together edges in breakpoint graph were indeed irregular
    (on end was incident to thevertex) infinity

    at the end a check is performed to make sure, that we've got the same number of genomes, from naming mapping,
    as we did from the lines, that contain gluing information

    Args:
        mgral_file: source text file with information about fusions. MGRAL file format is described in overview
        docstring

    Returns:
        tuple of two elements:
            1. dictionary of lists, where each genome (represented with one letter/digit alias) is mapped to a list
                of gluings, that were performed in it
            2. a dict of tuples, where each genome (represented by one letter/digit alias) is mapped to a a tuple,
                that contains 4 gene_name_prefix and genome verbose name
    """
    result_gluing_lists = defaultdict(list)
    result_genome_aliases = {}
    with open(mgral_file, "r") as source:
        for line in source:
            line = line.strip()
            if gnsregex.match(line):
                result = gnsregex.search(line)
                try:
                    one_letter_alias = result.groupdict()['ola']
                    genome_key_name = result.groupdict()['gkn']
                    genome_name = result.groupdict()['gn']
                    result_genome_aliases[one_letter_alias] = (genome_key_name, genome_name)
                except KeyError:
                    continue
            else:
                data = line.split("  ")
                for value in data:
                    if vrregex.match(value):
                        result = vrregex.search(value).groupdict()
                        try:
                            v1, v2, v3, v4 = result['vertex1'], result['vertex2'], result['vertex3'], result['vertex4']
                            ola = result['genome_alias']
                            if v1 and v2 or not v1 and not v2:
                                continue
                            if v3 and v4 or not v3 and not v4:
                                continue
                            gluing_pair = tuple(filter(lambda x: x is not None, (v1, v2, v3, v4)))
                            result_gluing_lists[ola].append(gluing_pair)
                        except KeyError:
                            continue
    assert set(x[0] for x in result_genome_aliases) == set(result_gluing_lists.keys())
    return result_genome_aliases, result_gluing_lists


def retrieve_genome_from_file(gff_file):
    """ reads gff file and returns a dict object, representing genome as a set of fragments

    retrieved data is split according to the sequence id, and information regarding "gene id", "start bp coordinate",
    "end bp coordinate", "strand' is retrieved.
    sorting of genes on each fragment is performed with respect to their base-pair coordinates
    tandem duplication is performed on each fragment, as we need to rewrite sequences of genes with just single
    instances of them for genome rearrangements purposes.

    Args:
        gff_file: text file in gff format. gff format is described in overview docstring

    Returns:
        a dict, where key -- name of fragment in genome (retrieved fro sequence id column in gff file), and value
        is a list of genes, that belong to this fragment.

    Raises:
        IndexError, as some lines might not go flawless through split process. As gff format is hugely used, any
            such inconsistency would terminate the process, as input data would be most likely corrupted and thus not
            reliable
    """
    result = defaultdict(list)
    with open(gff_file, "r") as source:
        for line in source:
            line = line.strip()
            data = line.split("\t")
            fragment = data[0]
            start, end = int(data[3]), int(data[4])
            strand = data[6]
            gene_id = data[10]
            result[fragment].append((gene_id, start, end, strand))
    for fragment in result:
        result[fragment] = sorted(result[fragment], key=lambda x: (x[1], x[2]))
        result[fragment] = remove_tandem_duplications(result[fragment])
    return result


def find_fragment_by_extremity(genome, gene_id, vertex):
    """ returns a name of a fragment and its direction, based on the supplied extremity

    iterates over fragments in genome and return the first fragment, that contains supplied gene on one of its two
    extremities.

    performs checks, to make sure that strand, gene is located on, corresponds to the end of a gene, according to the
    supplied value if breakpoint graph vertex

    Args:
        genome: a dict, where fragments names are mapped to a list of sorted genes
        gene_id: a gene_id, that suppose to represent an extremity of some of a fragment
        vertex: a gene_id representation in terms of breakpoint graph, must ends with "t" or "h"

    Returns:
    a tuple of two elements:
        1. fragment name, which contains supplied gene_id as its extremity
        2. direction of such fragment (or the end, on which such extremity is located:
            a. 1 if the fragment contains the extremity at its end
            b. -1 if the fragment contains the extremity at its start

    Raises:
        IndexError, if supplied genome, isn't represented as expected
        AssertionError: if orientation of a gene on the extremity of fragment doesn't correspond to its orientatiton
            in breakpoint graph
    """
    for fragment in genome:
        if genome[fragment][0][0] == gene_id:
            if len(genome[fragment]) == 1:
                if genome[fragment][0][3] == "+":
                    if vertex[-1] == "h":
                        return fragment, 1
                    else:
                        return fragment, -1
                else:
                    if vertex[-1] == "h":
                        return fragment, -1
                    else:
                        return fragment, 1
            if genome[fragment][0][3] == "+":
                assert vertex[-1] == "t"
            else:
                assert vertex[-1] == "h"
            return fragment, -1
        if genome[fragment][-1][0] == gene_id:
            if genome[fragment][-1][3] == "+":
                assert vertex[-1] == "h"
            else:
                assert vertex[-1] == "t"
            return fragment, 1


# TODO: implement algorithm
def connection_creation(f1, f2, d1, d2, storage):
    """ fills supplied storage with information, regarding connection between different fragments

    determines an extremity of each fragment (of supplied pair), based on directions, and append information
    to storage about connection, that particular fragment has in particular extremity

    Args:
        f1: name of first genome fragment
        f2: name of second genome fragment
        d1: direction (extremity) of fragment one, that was used in gluing operation
        d2: direction (extremity) of fragment two, that was used in gluing operation
        storage: data storage, that contains information about connections between different genome fragments

    Returns:
        Nothing, modifies storage object inplace
    """
    v1 = "h" if d1 == 1 else "t"
    v2 = "h" if d2 == -1 else "t"
    storage[f1][v1].append((f2, v2))
    storage[f2][v2].append((f1, v1))


# TODO: implement algorithm
def chain_construction(start_fragment, extremity, storage, visited):
    """ constructs a chain of connected fragments by pairwise connection data

    starting from supplied fragment, determines fragment, this one is connected to, and processes with it, appending
    resulted chain variable. If fragment has no further connection, resulted chain is returned.
    visited data storage object is updated, and double checked each time, algorithm proceeds, as only paths are possible
    in assembly, thus each fragment shall be visited exactly once

    Args:
        start_fragments: name of the fragment, that represents one of two ends of a future chain
        extremity: extremity of a fragment, that is connected to some other fragment
        storage: data storage object, that contains information about pairwise connection, between fragments
        visited: a data storage object, that contains information about fragments, that have been already traversed

    Returns:
        chain, a list of connected between each other fragments, with their orientation, in this chain
    """
    assert not visited[start_fragment]
    assert len(storage[start_fragment]) < 2

    visited[start_fragment] = True
    if extremity is None:
        return [(start_fragment, "+")]
    d = "+" if extremity == "h" else "-"
    result_chain = [(start_fragment, d)]
    extremity_to_go_to = "h" if extremity == "t" else "t"
    fragment = start_fragment
    while extremity_to_go_to in storage[fragment]:
        extremity_came_from = storage[fragment][extremity_to_go_to][1]
        extremity_to_go_to = "h" if extremity_came_from == "t" else "h"
        fragment = storage[fragment][extremity_to_go_to][0]
        visited[fragment] = True
        d = "+" if extremity_to_go_to == "h" else "-"
        result_chain.append((fragment, d))
    return result_chain


# TODO: implement algorithm
def get_assembly_fragments(genome, pairwise_gluing_info):
    """ rewrites genome fragment representation by assembling fragments together by pairwise gluing

    using information of pairwise gluing operations, algorithms performs construction of paths, of such glued fragments.
    once every fragment is visited through such chain creation, resulted value is return

    Args:
        genome: dictionary of lists of gene, where each key represents a genome fragment, and value stands for gene,
            that are located on this fragment
        pairwise_gluing_info:

    Returns:
        a set of lists, where each list stands for fragment chain (with length varying from 1 to ...)
    """
    connection_info_storage = defaultdict(lambda x: defaultdict(list))
    visited = {fragment_name: False for fragment_name in genome}
    for fragment_1, fragment_2, direction_1, direction_2 in pairwise_gluing_info:
        connection_creation(f1=fragment_1, f2=fragment_2, d1=direction_1, d2=direction_2,
                            storage=connection_info_storage)
    chains = set()
    ############################################################################################################
    #
    ############################################################################################################
    chain_started_fragments = filter(lambda fn: len(connection_info_storage[fn]) < 2, genome)
    for fragment_name in chain_started_fragments:
        if not visited[fragment_name]:
            ####################################################################################################
            #
            ####################################################################################################
            if len(connection_info_storage[fragment_name]) == 1:
                extremity = connection_info_storage[fragment_name].keys()[0]
            else:
                extremity = None
            chains.add(chain_construction(start_fragment=fragment_name,
                                          extremity=extremity,
                                          storage=connection_info_storage,
                                          visited=visited))
    return chains


def main(mgral_file, gene_mapping_file, gff_files):
    gene_number_mapping = retrieve_gene_number_mapping_from_file(gene_mapping_file=gene_mapping_file)
    genome_aliases, gluing_results = retrieve_glued_vertices_from_file(mgral_file=mgral_file)
    genomes = {}
    result = defaultdict(list)
    for gff_file in gff_files:
        genome_name = os.path.basename(gff_file).split(".")[0]
        genomes[genome_name] = retrieve_genome_from_file(gff_file)
    for genome_one_letter_alias in gluing_results:
        genome_name = genome_aliases[genome_one_letter_alias][1]
        genome_key = genome_aliases[genome_one_letter_alias][0]
        for v1, v2 in gluing_results[genome_one_letter_alias]:
            fragment1, d1 = find_fragment_by_extremity(genomes[genome_name],
                                                       gene_number_mapping[genome_key][int(v1[:-1])], v1)
            fragment2, d2 = find_fragment_by_extremity(genomes[genome_name],
                                                       gene_number_mapping[genome_key][int(v2[:-1])], v2)
            # assert fragment1 != fragment2
            result[genome_name].append((fragment1, fragment2, d1, d2))
    chained_result = {}
    for genome_name in genomes:
        chained_result[genome_name] = get_assembly_fragments(genome=genomes[genome_name],
                                                             pairwise_gluing_info=result[genome_name])
    # for genome in result:
    #     print(genome)
    #     gluing_statistics = result[genome]
    #     for f1, f2, d1, d2 in gluing_statistics:
    #         if f1 == f2:
    #             # print("\t", f1, "!!! (circularized)")
    #             continue
    #         if d1 == 1:
    #             if d2 == -1:
    #                 print("\t", f1, "(+) <==>", f2, "(+)")
    #             else:
    #                 print("\t", f1, "(+) <==>", f2, "(-)")
    #         else:
    #             if d2 == -1:
    #                 print("\t", f1, "(-) <==>", f2, "(+)")
    #             else:
    #                 print("\t", f2, "(+) <==>", f1, "(-)")


if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    main(cmd_args[0], cmd_args[1], cmd_args[2:])
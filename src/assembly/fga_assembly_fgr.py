#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from collections import defaultdict
import os
import re
import argparse

__author__ = 'Sergey Aganezov'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'

FRAGMENT_COLUMN = 0
START_COLUMN = 3
END_COLUMN = 4
STRAND_COLUMN = 6
GENE_ID_COLUMN = 10
SPLIT_CHAR = "\t"

# yeap, evil
# these regular expressions are used to parse MGRAL format file
vrregex = re.compile(
    "\((?:(?P<vertex1>\d+(?:h|t))|oo),(?:(?P<vertex2>\d+(?:h|t))|oo)\)x\((?:(?P<vertex3>\d+(?:h|t))|oo),(?:(?P<vertex4>\d+(?:h|t))|oo)\):\{(?P<genome_alias>\w)\}")
gnsregex = re.compile("(?P<ola>\w)=(?P<gkn>\w+)=(?P<gn>\w+)")


def retrieve_gene_id(line):
    """ Retrieves gene id from supplied string

    Uses script constants for determining which char use to split supplied line with and which item in split list to
    access

    Args:
        line: string of data

    Returns:
        gene id string

    Raises:
        IndexError: in case retrieving gene id from supplied line didn't go flawlessly
    """
    return line.split(SPLIT_CHAR)[GENE_ID_COLUMN]


def retrieve_start(line):
    """ Retrieves start base pair coordinate from supplied string

    Uses script constants for determining which char use to split supplied line with and which item in split list to
    access

    Args:
        line: string of data

    Returns:
        integer start bp coordinates

    Raises:
        IndexError: in case retrieving start coordinate from supplied line didn't go flawlessly
        ValueError: in case there was not an integer value in expected column in supplied string
    """
    return int(line.split(SPLIT_CHAR)[START_COLUMN])


def retrieve_end(line):
    """ Retrieves end base pair coordinate from supplied string

    Uses script constants for determining which char use to split supplied line with and which item in split list to
    access

    Args:
        line: string of data

    Returns:
        integer start bp coordinates

    Raises:
        IndexError: in case retrieving end coordinate from supplied line didn't go flawlessly
        ValueError: in case there was not an integer value in expected column in supplied string
    """
    return int(line.split(SPLIT_CHAR)[END_COLUMN])


def retrieve_fragment(line):
    """ Retrieves fragment name from supplied string

    Uses script constants for determining which char use to split supplied line with and which item in split list to
    access

    Args:
        line: string of data

    Returns:
        fragment name string

    Raises:
        IndexError: in case retrieving fragment name from supplied line didn't go flawlessly
    """
    return line.split(SPLIT_CHAR)[FRAGMENT_COLUMN]


def retrieve_strand(line):
    """ Retrieves DNA strand from supplied string

    Uses script constants for determining which char use to split supplied line with and which item in split list to
    access

    Args:
        line: string of data

    Returns:
        DNA strand (+|-) string

    Raises:
        IndexError: in case retrieving gene id from supplied line didn't go flawlessly
    """
    return line.split(SPLIT_CHAR)[STRAND_COLUMN]


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


def colapse_tandem_duplications(fragment_content):
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
    assert set(x[0] for x in result_genome_aliases) >= set(result_gluing_lists.keys())
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
            fragment = retrieve_fragment(line)
            start, end = retrieve_start(line), retrieve_end(line)
            strand = retrieve_strand(line)
            gene_id = retrieve_gene_id(line)
            result[fragment].append((gene_id, start, end, strand))
    for fragment in result:
        result[fragment] = sorted(result[fragment], key=lambda x: (x[1], x[2]))
        result[fragment] = colapse_tandem_duplications(result[fragment])
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
    v2 = "h" if d2 == 1 else "t"
    storage[f1][v1].append((f2, v2))
    storage[f2][v2].append((f1, v1))


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
    extremity_to_go_to = extremity
    fragment = start_fragment
    while extremity_to_go_to in storage[fragment]:
        extremity_came_from = storage[fragment][extremity_to_go_to][0][1]
        fragment = storage[fragment][extremity_to_go_to][0][0]
        extremity_to_go_to = "h" if extremity_came_from == "t" else "t"
        visited[fragment] = True
        d = "+" if extremity_to_go_to == "h" else "-"
        result_chain.append((fragment, d))
    return result_chain


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
    connection_info_storage = defaultdict(lambda: defaultdict(list))
    visited = {fragment_name: False for fragment_name in genome}
    for fragment_1, fragment_2, direction_1, direction_2 in pairwise_gluing_info:
        connection_creation(f1=fragment_1, f2=fragment_2, d1=direction_1, d2=direction_2,
                            storage=connection_info_storage)
    chains = []
    ############################################################################################################
    # at this point we try to take only those fragments, that have no more than one connection to some other fragment
    # in our sense these fragments are start/end of some chain (even knowing this chain can be of length 1)
    # chains are not double observed, as we maintain "visited" array to forbid that
    ############################################################################################################
    chain_started_fragments = filter(lambda fn: len(connection_info_storage[fn]) < 2, genome)
    for fragment_name in chain_started_fragments:
        if not visited[fragment_name]:
            ####################################################################################################
            # if we have indeed some connection to some other fragment, that means that this chains has at least
            # length 2, thus we have to which extremity on this tarting fragment we start chain retrieval
            ####################################################################################################
            if len(connection_info_storage[fragment_name]) == 1:
                extremity = list(connection_info_storage[fragment_name].keys())[0]
            else:
                ####################################################################################################
                # if there are no connections to other fragments, this means, that this fragment was not a part of
                # any gluing operation, and thus it represent a trivial chain of length 1
                ####################################################################################################
                extremity = None
            chains.append(chain_construction(start_fragment=fragment_name,
                                             extremity=extremity,
                                             storage=connection_info_storage,
                                             visited=visited))
    return chains


def main(mgral_file, gene_mapping_file, gff_files, options):
    gene_number_mapping = retrieve_gene_number_mapping_from_file(gene_mapping_file=gene_mapping_file)
    genome_aliases, gluing_results = retrieve_glued_vertices_from_file(mgral_file=mgral_file)
    genomes = {}
    result = defaultdict(list)
    for gff_file in gff_files:
        genome_name = os.path.basename(gff_file).split(".")[0].split("_")[0]
        genomes[genome_name] = retrieve_genome_from_file(gff_file)
    for genome_one_letter_alias in gluing_results:
        genome_name = genome_aliases[genome_one_letter_alias][1]
        genome_key = genome_aliases[genome_one_letter_alias][0]
        for v1, v2 in gluing_results[genome_one_letter_alias]:
            try:
                tmp_result = find_fragment_by_extremity(genomes[genome_name],
                                                        gene_number_mapping[genome_key][int(v1[:-1])], v1)
                if tmp_result is not None:  # None would be in the case, when supplied gene_id didn't correspond to any
                                            # outer-most synteny blocks in given genome
                    fragment1, d1 = tmp_result
                else:
                    continue
                tmp_result = find_fragment_by_extremity(genomes[genome_name],
                                                        gene_number_mapping[genome_key][int(v2[:-1])], v2)
                if tmp_result is not None:  # None would be in the case, when supplied gene_id didn't correspond to any
                                            # outer-most synteny blocks in given genome
                    fragment2, d2 = tmp_result
                else:
                    continue
            except AssertionError:  # an assertion error is raise, when supplied to find_fragment_by_extremity
                                    # gene id had incorrect orientation, in order for some extremity if it to take
                                    # place in a two-break operation
                continue
            result[genome_name].append((fragment1, fragment2, d1, d2))
    chained_result = {}
    for genome_name in genomes:
        if hasattr(options, "chains") and options.chains:
            chained_result[genome_name] = get_assembly_fragments(genome=genomes[genome_name],
                                                                 pairwise_gluing_info=result[genome_name])
        else:
            # if no chain option is supplied, each gluing will be represented a pair of fragments with respective
            # orientations
            visited = set()
            chained_result[genome_name] = []
            for value in result[genome_name]:
                fragment_1, fragment_2, d1, d2 = value
                # d1, d2 have values of 1 or -1, which has to be translated into + or - notation
                chained_result[genome_name].append([(fragment_1, "-" if d1 < 0 else "+"),
                                                    (fragment_2, "-" if d2 < 0 else "+")])
                visited.add(fragment_1)
                visited.add(fragment_2)
            for fragment in genomes[genome_name]:
                if fragment not in visited:  # if this fragment didn't take place in any gluing operation, it shall be
                                             #  added solely
                    chained_result[genome_name].append([(fragment, "+")])
        if hasattr(options, "genome_names") and options.genome_names:
            print(genome_name)
        for chain in chained_result[genome_name]:
            if hasattr(options, "glued") and options.glued:
                print("\t" + " <==> ".join("{f} ({d})".format(f=f, d=d) for f, d in chain))
            elif len(chain) > 1:
                print("\t" + " <==> ".join("{f} ({d})".format(f=f, d=d) for f, d in chain))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mgral_file", help="full path to mgra like log file with information about glued vertices")
    parser.add_argument("gene_num_mapping_file", help="full path to gene mapping file")
    parser.add_argument("-c", "--chains", action="store_true", default=False, dest="chains",
                        help="assembles glued pairs in glued chains")
    parser.add_argument("-g", "--glued", action="store_false", default=True, dest="glued",
                        help="permits output for assembled fragments only")
    parser.add_argument("gff_files", nargs="+", help="full path to gff formatted info files")
    parser.add_argument("--omit-genome-names", action="store_false", default=True, dest="genome_names",
                        help="omits printing of genome names in the output")
    args = parser.parse_args()
    main(mgral_file=args.mgral_file,
         gene_mapping_file=args.gene_num_mapping_file,
         gff_files=args.gff_files,
         options=args)
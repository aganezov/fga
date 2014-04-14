#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from collections import defaultdict
import sys

__author__ = 'Sergey Aganezov'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'

GENE_FAMILY_COLUMN = 1
ORGANISM_COLUMN = 4
GENE_ID_COLUMN = 3
SPLIT_CHAR = "\t"


def retrieve_gene_family(line):
    """ Retrieves gene family name from supplied string

    Uses script constants for determining which char use to split supplied line with and which item in split list to
    access

    Args:
        line: string of data

    Returns:
        gene family name string

    Raises:
        IndexError: in case retrieving gene family from supplied line didn't go flawlessly
    """
    return line.split(SPLIT_CHAR)[GENE_FAMILY_COLUMN]


def retrieve_organism(line):
    """ Retrieves organism name from supplied string

    Uses script constants for determining which char use to split supplied line with and which item in split list to
    access

    Args:
        line: string of data

    Returns:
        organism name string

    Raises:
        IndexError: in case retrieving organism name from supplied line didn't go flawlessly
    """
    return line.split(SPLIT_CHAR)[ORGANISM_COLUMN]


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


def get_set_of_organisms(data):
    """ Retrieves all organism names from supplied iterable

    For each entry in supplied iterable, retrieves organism and adds it to result set.
    Set is used to don't worry about duplications.

    Args:
        data: iterable, each entry in which is suitable for retrieve_organism function

    Returns:
        set of organisms name supplied iterable contains

    Raises:
        IndexError: as retrieve_organism might raise one
    """
    result = set()
    for entry in data:
        result.add(retrieve_organism(entry))
    return result


def get_unique_gene_families(data):
    """ Retrieves names of gene families that are present exactly once in each organism in supplied iterable

    First retrieves all existing organisms in supplied iterable
    Then for each gene family count how many times it occurs in each existing organism
    After that for each gene family we check that it occurs in ALL organisms and that it occurs in them no more than
    once

    Args:
        data: iterable, each entry in which is suitable for retrieve_organism, retrieve_gene_family functions

    Returns:
        set of gene family names, that occur exactly once in each organism in supplied iterable data

    Raises:
        IndexError: as get_set_of_organisms, retrieve_organism or retrieve_gene_family might raise one
    """
    gene_families = defaultdict(lambda: defaultdict(int))
    organisms = get_set_of_organisms(data)
    for entry in data:
        gene_families[retrieve_gene_family(entry)][retrieve_organism(entry)] += 1
    result = set()
    for gene_family in gene_families:
        if len(gene_families[gene_family].values()) == len(organisms) and \
                all(map(lambda x: x == 1, gene_families[gene_family].values())):
            result.add(gene_family)
    return result


def filter_gene_families(data, bad_gene_families=None, good_gene_families=None):
    """ scans through supplied iterable and check if retrieved gene family suites supplied filter sets

    function works as an iterator
    there are two supplied filter sets: good and bad once
    if good is present:
        function makes sure, that only those gene families, that present in good set are left
    if bad is present:
        function makes sure, that all gene families, that are present in bad set are filtered out

    bad set has priority over good one: if gene family name is present in both filter sets, it won't be left

    Args:
        data: iterable, each entry in which is suitable for retrieve_gene_family function
        bad_gene_families: set of gene family names, that must be filtered out
        good_gene_families: set of gene family names, that must be kept

    Returns:
        yields an entry in supplied iterable if it passed filtration

    Raises:
        IndexError: as retrieve_gene_family might raise one
    """
    for entry in data:
        gene_family = retrieve_gene_family(entry)
        good = True
        if bad_gene_families is not None and gene_family in bad_gene_families:
            good = False
        if good_gene_families is not None and gene_family not in good_gene_families:
            good = False
        if good:
            yield entry


def filter_organisms(data, good_organisms=None, bad_organisms=None):
    """ scans through supplied iterable and check if retrieved organism name suites supplied filter sets

    function works as an iterator
    there are two supplied filter sets: good and bad once
    if good is present:
        function makes sure, that only those organism names, that present in good set are left
    if bad is present:
        function makes sure, that all organism names, that are present in bad set are filtered out

    bad set has priority over good one: if organism name is present in both filter sets, it won't be left

    Args:
        data: iterable, each entry in which is suitable for retrieve_gene_family function
        bad_organisms: set of organism names, that must be filtered out
        good_organisms: set of organism names, that must be kept

    Returns:
        yields an entry in supplied iterable if it passed filtration

    Raises:
        IndexError: as retrieve_organism might raise one or if split of organism name won't go flawlessly
    """
    for entry in data:
        organism_name = retrieve_organism(entry)
        good = True
        if good_organisms is not None\
                and (organism_name not in good_organisms
                     and (len(organism_name.split(" ")) > 1 and organism_name.split(" ")[1] not in good_organisms)):
            good = False
        if bad_organisms is not None\
                and (organism_name in bad_organisms
                     or (len(organism_name.split(" ")) > 1 and organism_name.split(" ")[1] in bad_organisms)):
            good = False
        if good:
            yield entry


def main(mapping_file, bad_families_file=None, good_families_file=None,
         bad_organisms_file=None, good_organisms_file=None, settings=None):
    raw_data = mapping_file.readlines()
    ################################################################################################
    # some dirty tricks to work with both filtered and original ODB files, that contain first row with
    # column definition
    data = raw_data[1:] if raw_data[0].startswith("ODBMOZ") else raw_data
    data = list(filter(lambda y: y.split(SPLIT_CHAR)[1].startswith("MZ"), map(lambda x: x.strip(), data)))
    # we have to filter and leave only those rows, that contain legit info, which we assume for now, that every
    # row that starts with "MZ" string is legit. Will be changed in the future
    ################################################################################################
    if bad_families_file is not None or good_families_file is not None:
        bad_families, good_families = None, None
        if good_families_file is not None:
            good_families = set()
            for line in good_organisms_file:
                good_families.add(retrieve_gene_family(line.strip()))
        if bad_families_file is not None:
            bad_families = set()
            for line in bad_families_file:
                bad_families.add(line.strip())
        data = list(filter_gene_families(data=data, bad_gene_families=bad_families, good_gene_families=good_families))
    if good_organisms_file is not None or bad_organisms_file is not None:
        good_organisms, bad_organisms = None, None
        if good_organisms_file is not None:
            good_organisms = set()
            for line in good_organisms_file:
                good_organisms.add(line.strip())
        if bad_organisms_file is not None:
            bad_organisms = set()
            for line in bad_organisms_file:
                bad_organisms.add(line.strip())
        data = list(filter_organisms(data=data, good_organisms=good_organisms, bad_organisms=bad_organisms))
    if hasattr(settings, "unique") and settings.unique:
        good_families = get_unique_gene_families(data)
        data = list(filter_gene_families(data=data, good_gene_families=good_families))
    print("\n".join(data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("orthology_file", nargs="?", type=argparse.FileType("r"), default=sys.stdin,
                        help="full path to orthology files. Format can be seen in docs section,"
                             " in orth section. Or just standard input.")
    parser.add_argument("--bad-families-file", dest="bad_families_file", type=argparse.FileType("r"),
                        help="full file name with gene families ids to filter out",
                        default=None)
    parser.add_argument("--good-families-file", dest="good_families_file", type=argparse.FileType("r"),
                        help="full file name with gene families ids to keep during filtration",
                        default=None)
    parser.add_argument("--good-organisms-file", dest="good_organisms_file", type=argparse.FileType("r"),
                        help="full file name with organisms names to keep during filtration",
                        default=None)
    parser.add_argument("--bad-organisms-file", dest="bad_organisms_file", type=argparse.FileType("r"),
                        help="full file name with organisms names to filter out",
                        default=None)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-u", "--unique", action="store_true", default=False, help="Keeps only gene families, that are,"
                                                                                  " after filtration, present exactly"
                                                                                  " once in each organism")
    group.add_argument("-i", "--indels", help="doesn't work yet =(", action="store_true")
    group.add_argument("-d", "--duplications", help="doesn't work yet =(", action="store_true")
    args = parser.parse_args()
    main(args.orthology_file,
         bad_families_file=args.bad_families_file,
         good_families_file=args.good_families_file,
         bad_organisms_file=args.bad_organismss_file,
         good_organisms_file=args.good_organisms_file, settings=args)
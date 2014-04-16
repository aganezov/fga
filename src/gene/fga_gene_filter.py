#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from collections import defaultdict
import sys

__author__ = "Sergey Aganezov"
__email__ = "aganezov@gwu.edu"
__status__ = "develop"

FRAGMENT_COLUMN = 0
START_COLUMN = 3
END_COLUMN = 4
STRAND_COLUMN = 6
GENE_ID_COLUMN = 10
SPLIT_CHAR = "\t"


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


def filter_by_gene_ids(data, good_gene_ids=None, bad_gene_ids=None):
    """ Scans through supplied iterable and check if retrieved gene id suites supplied filter sets

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
        good = True
        gene_id = retrieve_gene_id(entry)
        if good_gene_ids is not None and gene_id not in good_gene_ids:
            good = False
        if bad_gene_ids is not None and gene_id in bad_gene_ids:
            good = False
        if good:
            yield entry


def get_genes_median_coordinates(data):
    """ For every gene id in supplied iterable data gets new coordinates, and returns obtained results

    retrieves all gene ids from supplied data
    for every gene id occurrence function computes median of its coordinates (start + end) / 2
    if there's more than one occurrence for some gene id, finish coordinates are computed as median of medians

    Args:
        data: iterable, each entry in which is suitable for retrieve_gene_id, retrieve_start, retrieve_end functions

    Returns:
        dict, where each key is a gene id, while value - new computed coordinate (median of median)

    Raises:
        IndexError: as retrieve_gene_id, retrieve_start, retrieve_end might raise one
        ValueError: as retrieve_start, retrieve_end might raise one
    """
    genes = defaultdict(list)
    result = {}
    for entry in data:
        gene_id = retrieve_gene_id(entry)
        start = retrieve_start(entry)
        end = retrieve_end(entry)
        genes[gene_id].append((start, end))
    for gene_id in genes:
        gene_cnt = len(genes[gene_id])
        gene_median_start = sum(start for start, end in genes[gene_id]) / gene_cnt
        gene_median_end = sum(end for start, end in genes[gene_id]) / gene_cnt
        result[gene_id] = (gene_median_start + gene_median_end) / 2
    return result


def filter_tandem_duplication(data):
    """ Removes repetitive fragments (by substituting them with only one member from them)

    function works as an iterator
    for each entry in supplied data, function determines if such entry has been already yielded, and if not, yields it

    Args:
        data: iterable, each entry in which is suitable for retrieve_gene_id function

    Returns:
        yields an entry in supplied iterable if its gene id value has not been yielded previously

    Raises:
        IndexError: as retrieve_gene_id might raise one
    """
    previous_gene_id = retrieve_gene_id(data[0])
    yield data[0]
    for entry in data[1:]:
        gene_id = retrieve_gene_id(entry)
        if gene_id != previous_gene_id:
            yield entry
            previous_gene_id = gene_id


# TODO: rewrite replacement part, as not only strings can be supplied
def rewrite_gene_coordinates(data, gene_id_coordinates_storage):
    """ Rewrites gene bp coordinates for every supplied entry in iterable data

    function works as an iterator
    for each entry in supplied data, function determines new coordinates for respective gene id and yields it

    Args:
        data: iterable, each entry in which is suitable for retrieve_gene_id, retrieve_start, retrieve_end function
        gene_id_coordinates_storage: dict, where key represents gene ids, and value new bp coordinates for respective
            gene ids

    Returns:
        yields an entry in which start and end coordinates are replaced with new ones, that were retrieved from
            supplied info storage

    Raises:
        IndexError: as retrieve_gene_id, retrieve_start, retrieve_end might raise one
    """
    for entry in data:
        start = retrieve_start(entry)
        end = retrieve_end(entry)
        gene_id = retrieve_gene_id(entry)
        entry = entry.replace("\t" + str(start) + "\t", "\t" + str(int(gene_id_coordinates_storage[gene_id])) + "\t")
        entry = entry.replace("\t" + str(end) + "\t", "\t" + str(int(gene_id_coordinates_storage[gene_id])) + "\t")
        yield entry


def retrieve_non_continuous_gene_ids(data):
    """ Scans among supplied data iterable and retrieves gene id names, that might be found in non contiguous sequences

    If gene ids appear as a non contiguous sequence in terms of bp coordinates, this function reports such gene id

    Args:
        data: iterable, that is suitable for retrieve_gene_id function

    Returns:
        set of gene ids, that were identified as being not just in contiguous sequences

    Raises:
        IndexError: as retrieve_gene_id might raise one
    """
    result = set()
    visited = {}
    previous_gene_id = retrieve_gene_id(data[0])
    visited[previous_gene_id] = True
    for entry in data[1:]:
        gene_id = retrieve_gene_id(entry)
        if gene_id != previous_gene_id:
            if gene_id in visited:
                result.add(gene_id)
            previous_gene_id = gene_id
            visited[gene_id] = True
    return result


def main(gff_file, good_gene_ids_file=None, bad_gene_ids_file=None, settings=None):
    data = gff_file.readlines()
    data = list(map(lambda x: x.strip(), data))
    if good_gene_ids_file is not None or bad_gene_ids_file is not None:
        good_gene_ids, bad_gene_ids = None, None
        if good_gene_ids_file is not None:
            good_gene_ids = set()
            for line in good_gene_ids_file:
                good_gene_ids.add(line.strip())
        if bad_gene_ids_file is not None:
            bad_gene_ids = set()
            for line in bad_gene_ids_file:
                bad_gene_ids.add(line.strip())
        data = list(filter_by_gene_ids(data, good_gene_ids=good_gene_ids, bad_gene_ids=bad_gene_ids))
    if hasattr(settings, "continuous") and settings.continuous:
        tmp_data = sorted(data,
                          key=lambda entry: (retrieve_fragment(entry), retrieve_start(entry), retrieve_end(entry)))
        non_continuous_gene_ids = retrieve_non_continuous_gene_ids(tmp_data)
        data = list(filter_by_gene_ids(data, bad_gene_ids=non_continuous_gene_ids))
    if hasattr(settings, "median") and settings.median:
        gene_median_coordinates = get_genes_median_coordinates(data)
        data = list(rewrite_gene_coordinates(data, gene_median_coordinates))
    if hasattr(settings, "sorted") and settings.sorted:
        data = sorted(data,
                      key=lambda entry: (retrieve_fragment(entry), retrieve_start(entry), retrieve_end(entry)))
    if hasattr(settings, "tandem_filtration") and settings.tandem_filtration:
        data = list(filter_tandem_duplication(data))
    print("\n".join(data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff_file", nargs="?", type=argparse.FileType("r"), default=sys.stdin,
                        help="full path to gff file. Format can be seen in docs section,"
                             " in gene section. Or just standard input.")
    parser.add_argument("--good-gene-ids-file", type=argparse.FileType("r"), help="full file name with gene ids to be"
                                                                                  " kept during filtration",
                        dest="good_gene_ids_file", default=None)
    parser.add_argument("--bad-gene-ids-file", type=argparse.FileType("r"), help="full file name with gene ids to"
                                                                                 " be filtered out during filtration",
                        dest="bad_gene_ids_file", default=None)
    parser.add_argument("--tandem-filtration", dest="tandem_filtration", action="store_true", default=False,
                        help="substitutes every tandem duplication of gene ids, with just one copy")
    parser.add_argument("--sorted", default=False, action="store_true",
                        dest="sorted", help="sorts all coding exons, that survived filtration, on respective fragments")
    parser.add_argument("-m", "--median", action="store_true", default=False,
                        help="rewrites each gene ids coordinates with median of median among all same "
                             "gene ids coordinates")
    parser.add_argument("-c", "--continuous", default=False, dest="continuous",
                        help="filters out all coding exons, that, if being sorted by bp coordinates,"
                             " don't form a contiguous sequence",
                        action="store_true")
    args = parser.parse_args()
    main(gff_file=args.gff_file, good_gene_ids_file=args.good_gene_ids_file, bad_gene_ids_file=args.bad_gene_ids_file,
         settings=args)

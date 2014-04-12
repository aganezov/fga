#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse
import sys
import itertools


__author__ = "Sergey Aganezov"
__email__ = "aganezov@gwu.edu"
__status__ = "develop"

GENE_FAMILY_COLUMN = 1
GENE_ID_COLUMN = 3
SPLIT_CHAR = "\t"

next_int = itertools.count()


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


def create_gene_family_number_mapping(data):
    """ Maps each gene family to some integer

    iterates over all gene families in supplied data iterable and maps them all to respective range of integers,
    starting with 0

    Args:
        data: iterable, where each entry suppose to contain gene family name and be processed by retrieve_gene_family
            function

    Returns:
        dist, where key stands for gene family name and value goes for integer it was mapped to

    Raises:
        IndexError: as retrieve_gene_family might raise one
    """
    gene_families = list(set(map(retrieve_gene_family, data)))
    gene_families = sorted(gene_families)
    return {gene_family: number for number, gene_family in enumerate(gene_families)}


def main(orth_file, settings):
    """ Performs and outputs gene id -> number mapping

    Reads supplied orthology mapping file and for each existing gene id in it outputs its mapping to some integer value,
     which was determined based on gene family particular gene id belongs to

    Args:
        orth_file: stream (file like object)

    Returns:
        nothing, writes to standard output

    Raises:
        IndexError, as retrieve_gene_family and retrieve_gene_id might raise one
    """
    if not settings.not_sorted:
        ####################################################################################
        # if we want to first sort gene families and then to map them to integers,
        # we have to use O(n) memory, and O(n*log(n)) time
        ####################################################################################
        data = orth_file.readlines()
        data =data[1:] if data[0].startswith("ODBMOZ") else data
        data = list(map(lambda l: l.strip(), data))
        gfnm = create_gene_family_number_mapping(data)
        for entry in data:
            gene_family = retrieve_gene_family(entry)
            gene_id = retrieve_gene_id(entry)
            print("{gene_id} {value}".format(gene_id=gene_id, value=gfnm[gene_family]))
    else:
        ####################################################################################
        # if we don't want to first sort gene families and then to map them to integers,
        # we have to use O(1) memory and O(n) time
        ####################################################################################
        gfnm = {}
        for line in orth_file:
            if line.startswith("ODBMOZ"):
                continue
            gene_family = retrieve_gene_family(line)
            gene_id = retrieve_gene_id(line)
            if gene_family not in gfnm:
                gfnm[gene_family] = next(next_int)
            print("{gene_id} {value}".format(gene_id=gene_id, value=gfnm[gene_family]))


if __name__ == "__main__":
    python_version = sys.version_info
    if python_version.major < 3:
        print("This script is suitable for only python 3.3+", file=sys.stderr)
        print("Your python version is", ",".join(map(str, python_version)), file=sys.stderr)
        exit("-1")
    if python_version.major == 3 and python_version.minor < 3:
        print("This script is suitable for only python 3.3+", file=sys.stderr)
        print("Your python version is", "-".join(map(str, python_version)), file=sys.stderr)
        exit("-1")
    parser = argparse.ArgumentParser()
    parser.add_argument("orth_file", nargs="?", type=argparse.FileType("r"), default=sys.stdin,
                        help="full path to orthology file. Format can be seen in docs section,"
                             " in orth section. Or just standard input.")
    parser.add_argument("--not-sorted", dest="not_sorted", action="store_false", default=True,
                        help="first sorts all gene families and then maps them to integers, requires O(n) memory,"
                             "where n equals to the number of line in source")
    args = parser.parse_args()
    main(orth_file=args.orth_file, settings=args)

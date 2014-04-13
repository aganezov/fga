#!/usr/local/bin/python3.4
# -*- coding: utf-8 -*-
import argparse
from collections import defaultdict
import os
import sys

__author__ = "Sergey Aganezov"
__email__ = "aganezov@gwu.edu"
__status__ = "develop"

FRAGMENT_COLUMN = 0
START_COLUMN = 3
END_COLUMN = 4
STRAND_COLUMN = 6
GFF_GENE_ID_COLUMN = 10
GFF_SPLIT_CHAR = "\t"
GNM_SPLIT_CHAR = " "
GNM_GENE_ID_COLUMN = 0
GNM_MAPPING_VALUE_NUMBER = 1


def retrieve_gene_id_gff(line):
    """ Retrieves gene id from supplied string using knowledge of gff files format

    Uses script constants for determining which char use to split supplied line with and which item in split list to
    access

    Args:
        line: string of data

    Returns:
        gene id string

    Raises:
        IndexError: in case retrieving gene id from supplied line didn't go flawlessly
    """
    return line.split(GFF_SPLIT_CHAR)[GFF_GENE_ID_COLUMN]


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
    return int(line.split(GFF_SPLIT_CHAR)[START_COLUMN])


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
    return int(line.split(GFF_SPLIT_CHAR)[END_COLUMN])


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
    return line.split(GFF_SPLIT_CHAR)[FRAGMENT_COLUMN]


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
    return line.split(GFF_SPLIT_CHAR)[STRAND_COLUMN]


def retrieve_gene_id_gnm(line):
    """ Retrieves gene id from supplied string using knowledge of gene_number_mapping files format

    Uses script constants for determining which char use to split supplied line with and which item in split list to
    access

    Args:
        line: string of data

    Returns:
        gene id string

    Raises:
        IndexError: in case retrieving gene id from supplied line didn't go flawlessly
    """
    return line.split(GNM_SPLIT_CHAR)[GNM_GENE_ID_COLUMN]


def retrieve_value(line):
    """ Retrieves mapping value for gene id from supplied string using knowledge of gene_number_mapping files format

    Uses script constants for determining which char use to split supplied line with and which item in split list to
    access

    Args:
        line: string of data

    Returns:
        value for gene id in respective mapping line string

    Raises:
        IndexError: in case retrieving value for gene id from supplied line didn't go flawlessly
    """
    return line.split(GNM_SPLIT_CHAR)[GNM_MAPPING_VALUE_NUMBER]


def get_gnm(data):
    """ Retrieves overall gene ids -> number mapping from supplied iterable

    Args:
        data: iterable, each entry in which is suitable for retrieve_gene_id_gnm, retrieve_value functions

    Returns:
        dict:
            keys stand for gene ids
            values stand  values, that respective gene ids are mapped in supplied data iterable

    Raises:
        IndexError: as retrieve_gene_id_gnm or retrieve_value might raise one
    """
    return {retrieve_gene_id_gnm(line): retrieve_value(line) for line in data}


def retrieve_genome_name_from_file(gff_file):
    """ Returns genome_name from supplied file name

    usually genome file basename might contain _ and followed by some additional info,, thus only genome name shall be
    retrieved

    Args:
        fgg_file: full name for gff file

    Returns:
        genome name retrieved from supplied file name
    """
    return os.path.basename(gff_file).split(".")[0].split("_")[0]


def get_mgra_representation_of_a_gene(gene, gnm):
    """ Returns a string representation as oriented value for supplied gene

    mgra uses simple integer representation (without "+" char) for positively oriented blocks,
    while requires "-" char in front of synteny blocks from negative strand

    Args:
        gene: iterable value
            first element used for mapping to string representation in gnm
            fourth element used for determining orientation of respective block
        gnm: dict, using for mapping supplied gene value for its future string representation

    Returns:
        string representation for supplied gene value, taking orientation into consideration

    Raises:
        IndexError: as accessing supplied value by index might be costly =(
    """
    result = ""
    if gene[3] == "-":
        result += "-"
    result += str(gnm[gene[0]])
    return result


def read_genome_from_iterable(data):
    """ Reads a genome information from supplied iterable

    Genome is represented as a set of fragments, each of which is represented as a sequence of genes
    fragments are represented as keys in result dict
    sequences are represented with lists in resulted dict, which are mapped to fragment names, they belong to

    Args:
        data: iterable,
            each entry of which contains information about one gene in genome
            each entry of which must be suitable for retrieve_fragment, retrieve_gene_id_gff, retrieve_start,
                retrieve_end, retrieve_strand functions
    Returns:
        defaultdict of lists, where keys stand for fragment names in genome, while values represent sequences of genes
            on each fragments

    Raises:
        IndexError: as retrieve_fragment, retrieve_gene_id_gff, retrieve_start, retrieve_end or retrieve_strand might
            raise one
    """
    result = defaultdict(list)
    for entry in data:
        fragment = retrieve_fragment(entry)
        gene_id = retrieve_gene_id_gff(entry)
        start = retrieve_start(entry)
        end = retrieve_end(entry)
        strand = retrieve_strand(entry)
        result[fragment].append((gene_id, start, end, strand))
    return result


def main(gnm_file, gff_files):
    data = gnm_file.readlines()
    data = list(map(lambda l: l.strip(), data))
    gnm = get_gnm(data)
    for gff_file in gff_files:
        genome_name = retrieve_genome_name_from_file(gff_file)
        with open(gff_file, "r") as source:
            print(">" + str(genome_name))
            data = source.readlines()
            data = list(map(lambda l: l.strip(), data))
            genome = read_genome_from_iterable(data)
            for fragment in sorted(genome):
                genes = genome[fragment]
                genes = sorted(genes, key=lambda g: (g[1], g[2]))
                genes = list(map(lambda gene: get_mgra_representation_of_a_gene(gene, gnm), genes))
                print(" ".join(genes), "$")
        print()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gene_number_mapping_file", type=argparse.FileType("r"),
                        help="full path to gene number mapping file. Format can be seen in docs section,"
                             "in gene section.")
    parser.add_argument("gff_files", nargs="+",
                        help="full path to gff files for different genomes")
    args = parser.parse_args()
    main(args.gene_number_mapping_file, gff_files=args.gff_files)
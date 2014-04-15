# Documentation

This is a documentation section for **fga** project software, which will help you navigate around source code,
understand basic principals of what which script does, and hopefully answer all you questions.

## Purpose

The **FGA** project was undertaken to study evolutionary genomics in several mosquito genomes form South Africa.

As the data had to be prepared and analyzed, several software scripts for those purposes were created. All scripts
were written to be used as command line stand alone tools.

In general, all software is divided into following categories:

1. **orth** - this set of scripts deals with data, that contains information about orthology relationship among genomes
    (gene families). More on particular script usage and data format in [orthology section](#Orthology)
2. **gene** - this set of scripts deals with data, that contains information about coding exons and fragments in each
    genome. More on particular script usage and data format in [gene section](#Gene)
3. **assembly** - this set of scripts deals with data, that contains information about assembled statistics, among genomes
    More on particular script usage and data format in [assembly section](assembly)

Each set is contained in a respectively named python package. All scripts are written to be invoked as standalone command
line tools. All scripts write to standard output, thus can be used as unix pipes sources. Some scripts can read from standard
input only, if needed, thus can be used in any position in unix pipes.

All scripts are written with support of python3.3+. Each script contains ``#! /usr/bin/env python3`` directive, thus can be
invoked without using ``python3 script_name.py`` prefix.

All scripts are named as follows:

    <project_name>_<category_name>_<script_name>.py

As each script was meant to be standalone and not dependent on any other package, library (except standard), etc., code
contains some duplication in terms of support functions.


# Orthology

This section covers all scripts that can be found in ``orth`` python package in root ``src`` package. This scripts are
designed to process orthologous data, that is presented in tab-separated text format. All scripts assumes some things about
supplied data, but all those bottlenecks can be easily generalized.

Assumptions:


* all supplied data is separated by ``\t` symbols
* second column in supplied data contains information regarding orthologous families names
* forth column in supplied data contains information regarding particular genes, each gene family contains
* fifth column in cupplied data contains information regarding organisms name

Thus one row in supplied data read as follows:

    gene family ... gene id ... organism name

and this means, that organism ``organism name`` contains a gene ``gene id``, which belongs to gene family ``gene family``

As sometimes this data, while been contained in files, contains first-row-column definition, scripts skips first row in
 analysis.

## [fga_orth_filter.py][1]

This script is designed for filtration purposes, as most of the times not all existing raw data can be suited for experiment
Using this script one can apply filters to ``gene family`` and ``organism name`` columns in supplied data. Scripts outputs
rows from supplied data as is, is they suite filtration conditions.

#### Usage

    fga_orth_filter.py ODB_file

    cat ODB_file | fga_orth_filter.py

Such call will output supplied file as is (assuming all rows were suitable for processing), as no options were used for filtration.
This script allows one to supply text files with **good** and **bad** values for both filtration columns.

    fga_orth_filter.py ODB_file -u

This call will limit all supplied data only to rows that contain gene families, that were present exactly once in all
organisms that could have been found in supplied data.

    fga_orth_filter.py -u --good-organisms-file g_organisms_file.txt

This call will limit all supplied data only to rows that contain gene families, that were present exactly once in all
organisms that could have been found in supplied ``g_organisms_file.txt`` file.

As ``-h, --help`` option states:

    usage: fga_orth_filter.py [-h] [--bad-families-file BAD_FAMILIES_FILE]
                          [--good-families-file GOOD_FAMILIES_FILE]
                          [--good-organisms-file GOOD_ORGANISMS_FILE]
                          [--bad-organisms-file BAD_ORGANISMS_FILE]
                          [-u | -i | -d]
                          [orthology_file]

    positional arguments:
      orthology_file        full path to orthology files. Format can be seen in
                            docs section, in orth section. Or just standard input.

    optional arguments:
      -h, --help            show this help message and exit
      --bad-families-file BAD_FAMILIES_FILE
                            full file name with gene families ids to filter out
      --good-families-file GOOD_FAMILIES_FILE
                            full file name with gene families ids to keep during
                            filtration
      --good-organisms-file GOOD_ORGANISMS_FILE
                            full file name with organisms names to keep during
                            filtration
      --bad-organisms-file BAD_ORGANISMS_FILE
                            full file name with organisms names to filter out
      -u, --unique          Keeps only gene families, that are, after filtration,
                            present exactly once in each organism
      -i, --indels          doesn't work yet =(
      -d, --duplications    doesn't work yet =(

for this script to work there must be either a standard input stream, or full path to file, that contains data for
 filtration.

It's important to note, that **bad** filtration sets always dominate good ones: if there's a gene family in both good
and bad gene family sets, it will be discarded during the filtration process.

This script allows one to limit gene families to different sets, which might suite for different genome rearrangement
 analysis scenarios:

1. Analysis of genomes, that undertook evolutionary scenarios, that contain only genome rearrangement events.
 Suitable output for this purpose can be achieved by using **``-u``** option, when invoking the filtration script.
2. Analysis of genomes, that undertook evolutionary scenarios, that contain genome rearrangements, insertions and
 deletions of unique genomic content. ***Not yet implemented.***
3. Analysis of genomes, that undertook evolutionary scenarios, that contain genome rearrangements, insertions,
 deletions and duplication events. ***Not yet implemented.***

[1]:https://github.com/sergey-aganezov-jr/fga/blob/master/src/orth/fga_orth_filter.py

## [fga_orth_gnm.py][2]

This script is designed for creating so called gene-number mapping. As one gene family can contain multiple different
gene coding sequences in different genomes, all those genes, for the genome rearrangements shall be represented as same
 instances of homologous gene families. Thus each gene family gets a particular value (integer), and then all gene ids,
  that belong to this gene family would be assigned respective value and output.

#### Usage

    fga_orth_gnm.py ODB_file

    cat ODB_file | fga_orth_gnm.py

Such call will take all supplied data in, filter all gene families id, that could have been found in supplied data,
and then process supplied data in the **supplied** order, assigning each gene id to respective value, gene family, it's
a part of, was previously assigned to.

    fga_orth_filter.py ODB_file --non-sorted

Such call will stream supplied data as is. The ``--non-sorted`` option was disabled by default, as it comes more often
than not, to have consistency in terms of that mapping in regards of supplied data.

As ``-h/--help`` option states:

    usage: fga_orth_gnm.py [-h] [--non-sorted] [orth_file]

    positional arguments:
      orth_file     full path to orthology file. Format can be seen in docs
                    section, in orth section. Or just standard input.

    optional arguments:
      -h, --help    show this help message and exit
      --non-sorted  streams all rows from supplied stream/file as is, whithout
                    first filtering them by gene family name first (disabled by
                    default)

for this script there must be either an input stream, or a full path to file, that contains data for creating gene-number
mapping.

[2]:https://github.com/sergey-aganezov-jr/fga/blob/master/src/orth/fga_orth_gnm.py

## [fga_orth_check.py][3]

This script is designed for simple pairwise check of different orthology sets. In other terms, it check if same
genes from different organisms grouped in same families across different data files.

The main idea is as follows: for each gene family in first supplied file, determine which genes from which organisms,
 are mapped to this family. Assign same value to gene family and all genes that belong to it. For the second file do
 the same procedure, derive data in different chunks of gene families. For each gene family get all gene that belong
 to it in second data set. Get a value for this gene family as the most encountered in terms of what genes in this
 family were mapped in first file. Once the value is derived, check if there some gene in observed family, that were
 mapped differently. If there's at least one so gene - report the gene family, as the one having a miss-map.

#### Usage

    fga_orth_check.py ODB_file_1 ODB_file_2

Such call will take all supplied data in, determine gene family <-> gene id mapping based on first file, and then check
gene families in the second file. The output of this script is value from second column (gene families names) in file
``ODB_file_2`` that have at least one gene id miss mapped based on file one. So if one want to perform same analysis,
 based on different orhtology datasets, and wants to be sure in that initial data is consistent among those datasets,
 they should perform the following command sequence to collect all unreliable gene families:

    fga_orth_check.py ODB_file_1 ODB_file_2 > bad_gene_families

    fga_orth_check.py ODB_file_2 ODB_file_1 >> bad_gene_families

these calls will collect all gene families names, that contain anything inconsistent in them. Using
``fga_orth_filter.py`` after this with ``--bad-families-file`` option pointing to ``bad_gene_families`` file on both
 those ODB files would produce consistent with each other orthology mapping files.

As ``-h/--help`` option states:

    usage: fga_orth_check.py [-h] first_orthology_file second_orthology_mapping

    positional arguments:
      first_orthology_file  Full path to first orthology mapping file
      second_orthology_mapping
                            Full path to second orthology mapping file to check
                            against first

    optional arguments:
      -h, --help            show this help message and exit

this script expects exactly two arguments, which suppose to provide full paths to orthology mapping files, where second
 file has to be checked against the first one.

[3]:https://github.com/sergey-aganezov-jr/fga/blob/master/src/orth/fga_orth_check.py


# Gene

This section covers all scripts that can be found in ``gene`` python package in root ``src`` package. This scripts are
designed to process gff formatted data, that is presented in gff format. All scripts assumes that gff formatted files
have only correct info. More of gff format can be read [here](http://www.ensembl.org/info/website/upload/gff.html).
The only assumption, is that on the 11th column we'll have a **gene id** information, which is used wildly in these
 scripts.

When we talk about **gene <-> number mapping**  files, we assume the following format:

    gene_id integer_value

We assume two space-separated words, of which second consists of only digit (as it will be casted to integer value)


## [fga_gene_filter.py][4]

This script is designed for simple filtration of gff formatted files, which in this study were representing information
about particular gene ids location in each genome. There are several filtration available with this script:

* **good** / **bad** gene ids filter: as files with those values can be supplied for this script, one can filter out all
unwanted gene ids from supplied gff formatted data. It is important to notice, that **bad** filtration set will always
 dominate a good one.

As same gene can be represented as a sequence of coding exons, and, for the purpose of genome rearrangement perspective,
 we'd like to represent each gene with a single coding exon. There are several options available for rewriting coding
  exons:

1. median coordinates rewriting: for each gene we determine call coding exons. After that, for each coding exon
 we determine a median base pair coordinate. And the resulting coordinate for particular gene is computed as median of
  previously computed medians. Results are written into start and end coordinates of respective coding exon sequences.

2. sorting: on all coding fragments (first column sequences ids) script performs sorting of coding exons, using their
start and end coordinates.

3. tandem duplication elimination: if, there are several coding sequences are following each other, they'll be
represented with a single instance in the resulted data

in terms of order, filtration and rewriting is performed in the following order:

1. gene id filtration, using **bad** and **good** gene id sets
2. median coordinates rewriting
3. bp coordinates sorting
4. tandem duplication filtration

#### Usage

    fga_gene_filter.py gff_file

    cat gff_file | fga_gene_filter.py

Such call will simply output supplied data in gff format.

As ``-h/--help`` options states:

    usage: fga_gene_filter.py [-h] [--good-gene-ids-file GOOD_GENE_IDS_FILE]
                          [--bad-gene-ids-file BAD_GENE_IDS_FILE]
                          [--no-tandem-filtration] [--non-sorted] [-m] [-c]
                          [gff_file]

    positional arguments:
      gff_file              full path to gff file. Format can be seen in docs
                            section, in gene section. Or just standard input.

    optional arguments:
      -h, --help            show this help message and exit
      --good-gene-ids-file GOOD_GENE_IDS_FILE
                            full file name with gene ids to be kept during
                            filtration
      --bad-gene-ids-file BAD_GENE_IDS_FILE
                            full file name with gene ids to be filtered out during
                            filtration
      --no-tandem-filtration
                            stops substituting of every tandem duplication of gene
                            ids, with just one copy
      --not-sorted          prevents sorting of coding exons on respective
                            fragments
      -m, --median          rewrites each gene ids coordinates with median of
                            median among all same gene ids coordinates
      -c, --continuous      doesn't work yet =(

the script expect a full path to gff formatted file, or standard input. Byt the options flags, one can see, that bp
coordinates sorting is enabled by default, tandem filtration is enabled by default. While median coordinates rewriting
has to be invoked explicitly.

[4]:https://github.com/sergey-aganezov-jr/fga/blob/master/src/gene/fga_gene_filter.py
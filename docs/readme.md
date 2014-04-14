# Documentation

This is a documentation section for **fga** project software, which will help you navigate around source code,
understand basic principals of what which script does, and hopefully answer all you questions.

## Purpose

The **FGA** project was undertaken to study evolutionary genomics in several mosquito genomes form South Africa.

As the data had to be prepared and analyzed, several software scripts for those purposes were created. All scripts
were written to be used as command line stand alone tools.

In general, all software is divided into following categories:

1. **orth** - this set of scripts deals with data, that contains information about orthology relationship among genomes
    (gene families). More on particular script usage and data format in [orthology section](#orth)
2. **gene** - this set of scripts deals with data, that contains information about coding exons and fragments in each
    genome. More on particular script usage and data format in [gene section](#gene)
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


# <a id="orth"></a>Orthology

This section covers all scripts that can be found in ``orth`` python package in root ``src`` package. This scripts are
designed to process orthologous data, that is presented in tab-separated text format. All scripts assumes some thigngs about
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

    usage: fga_orth_gnm.py [-h] [--non-sorted] [orth_file]

    positional arguments:
      orth_file     full path to orthology file. Format can be seen in docs
                    section, in orth section. Or just standard input.

    optional arguments:
      -h, --help    show this help message and exit
      --non-sorted  streams all rows from supplied stream/file as is, whithout
                    first filtering them by gene family name first (disabled by
                    default)

[2]:https://github.com/sergey-aganezov-jr/fga/blob/master/src/orth/fga_orth_gnm.py
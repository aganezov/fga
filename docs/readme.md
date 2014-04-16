# Documentation

This is a documentation section for **fga** project software, which will help you navigate around source code,
understand basic principals of what which script does, and hopefully answer all you questions.

## Purpose

The **FGA** project was undertaken to study evolutionary genomics in several mosquito genomes form South Africa.

As the data had to be prepared and analyzed, several software scripts for those purposes were created. All scripts
were written to be used as command line standalone tools.

In general, all software is divided into following categories:

1. **orth** - this set of scripts deals with data, that contains information about orthology relationship among genomes
    (gene families). More on particular script usage and data format in [orthology section](https://github.com/sergey-aganezov-jr/fga/tree/master/docs#orthology)
2. **gene** - this set of scripts deals with data, that contains information about coding exons and fragments in each
    genome. More on particular script usage and data format in [gene section](https://github.com/sergey-aganezov-jr/fga/tree/master/docs#gene)
3. **assembly** - this set of scripts deals with data, that contains information about assembled statistics, among genomes
    More on particular script usage and data format in [assembly section](#assembly)

Each set is located in a respectively named python package in ``src`` package. All scripts are written to be invoked as
 standalone command line tools. All scripts write to standard output, thus can be used as unix pipes sources.
  Some scripts can read from standard input only, if needed, thus can be used in any position in unix pipes.

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


* all supplied data is separated by ``\t`` symbols
* second column in supplied data contains information regarding orthologous gene families names
* forth column in supplied data contains information regarding particular genes (gene ids, coding exons), each gene family contains
* fifth column in supplied data contains information regarding organisms name

Thus one row in supplied data read as follows:

    ... gene family ... gene id ... organism name ...

and this means, that organism ``organism name`` contains a gene ``gene id``, which belongs to gene family ``gene family``

As sometimes this data, while been contained in files, contains column definition in the first row, scripts skips
first row in analysis.

## [fga_orth_filter.py][1]

This script is designed for filtration purposes, as most of the times not all existing raw data can be suited for experiment.
Using this script one can apply filters to ``gene family`` and ``organism name`` columns in supplied data. Scripts outputs
rows from supplied data as is, if they suite filtration conditions.

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

1. Analysis of genomes, which undertook evolutionary scenarios, that contain only genome rearrangement events.
 Suitable output for this purpose can be achieved by using **``-u/--unique``** option, when invoking the filtration script.
2. Analysis of genomes, which undertook evolutionary scenarios, that contain genome rearrangements, insertions and
 deletions of unique genomic content. ***Not yet implemented.***
3. Analysis of genomes, which undertook evolutionary scenarios, that contain genome rearrangements, insertions,
 deletions and duplication events. ***Not yet implemented.***

 As one could have seen, from ``-h/--help`` option output, these scenarios are exclusive, thus only one option can be
 legitimately used at a time.

[1]:https://github.com/sergey-aganezov-jr/fga/blob/master/src/orth/fga_orth_filter.py

## [fga_orth_gnm.py][2]

This script is designed for creating so called gene-number mapping. As one gene family can contain multiple different
gene ids in different genomes, all those genes ids for the genome rearrangements purposes shall be represented as same
 instances of respective gene families. Thus each gene family gets a particular value (integer), and then all gene ids,
  that belong to this gene family would be assigned to a respective value.

#### Usage

    fga_orth_gnm.py ODB_file

    cat ODB_file | fga_orth_gnm.py

Such call will take all supplied data in, sort all gene families, that could have been found in supplied data,
and then process supplied data in the **supplied** order, assigning each gene id to respective value, gene family, the gene id is
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
genes from different organisms were grouped in same gene families across different data files.

The main idea is as follows: for each gene family in first supplied file, determine which genes from which organisms,
 are mapped to this family. Assign same value to gene family and all genes that belong to it. For the second file,
  derive data in different chunks of gene families. For each gene family get all gene that belong
 to it in second data set. Get a value for this gene family as the most encountered in terms of what genes in this
 family were mapped in first file. Once the value is derived, check if there some gene in observed gene family in second dataset,
  that were mapped differently in the first dataset. If there's at least one such miss-mapped gene - report the gene family, as the one having a miss-map.

#### Usage

    fga_orth_check.py ODB_file_1 ODB_file_2

Such call will take all supplied data in, determine gene family <-> gene id mapping based on first file, and then check
gene families in the second file. The output of this script will be values from second column (gene families names) in file
``ODB_file_2`` that have at least one gene id miss mapped based on file ``ODB_file_1``. So if one want to perform same analysis,
 based on different orhtology datasets, and wants to be sure that initial data is consistent among those datasets,
 they should perform the following command sequence to collect all unreliable gene families:

    fga_orth_check.py ODB_file_1 ODB_file_2 > bad_gene_families

    fga_orth_check.py ODB_file_2 ODB_file_1 >> bad_gene_families

these calls will collect all gene families names, that contain anything inconsistent in them. Usage of
``fga_orth_filter.py`` script after this, with ``--bad-families-file`` option pointing to ``bad_gene_families`` file, on both
 those ODB files, would produce consistent with each other orthology mapping files.

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
designed to process gff formatted data. All scripts assumes that gff formatted files
have only correctly formatted info. More of gff format can be read [here](http://www.ensembl.org/info/website/upload/gff.html).
The only assumption, is that on the 11th column, in gff formatted files, we'll have a **gene id** information, which is used wildly in these
 scripts and in ``orth`` scripts.

When we talk about **gene <-> number mapping**  files, we assume the following format:

    gene_id integer_value

We assume two space-separated words, of which second consists of only digit (as it will be casted to integer value)


## [fga_gene_filter.py][4]

This script is designed for simple filtration of gff formatted files, which in this study were representing information
about particular gene ids location in each genome. There are several filtrations available with this script:

* **good** / **bad** gene ids filter: as files with those values can be supplied for this script, one can filter out all
unwanted gene ids from supplied gff formatted data. It is important to notice, that **bad** filtration set will always
 dominate a good one.

* **continuous** gene id filter: this option allows one to examine all gene ids, as if they were sorted by their bp
 coordinates, and then to check if there are any gene ids, that might appear in non contiguous sequences. Such
 gene ids will be removed from further filtration.

As some gene can be represented as a sequence gene ids, and, from the genome rearrangement perspective,
 we'd like to represent each gene with a single gene. There are several options available for rewriting gene ids:

1. **median** coordinates rewriting: for each gene we determine call coding exons. After that, for each coding exon
     we determine a median base pair coordinate. And the resulting coordinate for particular gene is computed as median of
      previously computed medians. Results are written into start and end coordinates of respective coding exon sequences.


        Example:
        assume one has 3 identical gene ids AAA111 in some genome on same fragment with coordinates:

        1. start = 0, end = 5
        2. start = 20, end = 25
        3. start = 40, end = 45

        this option will rewrite for each of those gene ids their start and end coordinates as 22.5, as it equals to

        ((0 + 5) / 2 + (20 + 25) / 2 + (40 + 45) / 2 ) / 3

        co called median of median

2. **tandem filtration**: if, there are several gene ids, that are following each other in terms of their bp coordinates,
  they'll be represented with a single instance in the resulted data


        Example:
        assume one has 3 identical gene ids AAA111 in some genome on same fragment that form a contiguous sequence (as
         there are no other gene ids, that, based on their coordinates would be in between AAA111 ids after sorting)

        ... AAA111 AAA111 AAA111 ...

        this option will rewrite such gene ids as

        ... AAA111 ...

        leaving only one instance from each such contiguous sequence of same gene ids

in terms of order, filtration and rewriting are performed in the following order:

1. gene id filtration, using **bad** and **good** gene id sets
2. continuous filtration
3. median coordinates rewriting
4. bp coordinates sorting
5. tandem duplication filtration

#### Usage

    fga_gene_filter.py gff_file

    cat gff_file | fga_gene_filter.py

Such call will simply output supplied data in gff format.

    fga_gene_filter.py gff_file --bad-gene-ids-file bed_gene_ids

Such call will output supplied data, skipping all rows, where 11th column contains values, that could be found in ``bed_gene_ids`` file

    fga_gene_filter.py gff_file -c

    fga_gene_filter.py gff_file --continuous

Such call will take all supplied data in, internally sort it, determine which gene ids could have been found in non contiguous sequences,
and then output supplied data in **supplied** order, skipping those rows, where those non-contiguous gene ids were specified in 11th column

    fga_gene_filter.py gff_file -c -m

    fga_gene_filter.py gff_file -c | fga_gene_filter.py -m

Such call will first filter out all rows which contain non-contiguous gene ids, and then will rewrite coordinates for all
gene ids in all rows, with newly computed median of median values.

    fga_gene_filter.py gff_file -m --sorted --tandem-filtration

    fga_gene_filter.py gff_file -m | fga_gene_filter.py gff_file --sorted | fga_gene_filter.py gff_file --tandem-filtration

Such call will for all gene ids rewrite their coordinates, then sort them, using newly computed coordinates and then collapse
contiguous sequences of same gene ids into single instances of those repetitive gene ids.

**IMPORTANT** without specifying ``--sorted`` key, the result won't be correct, as while coordinates will be rewritten for
 all gene ids, tandem filtration will be performed based on **supplied** order.

    fga_gene_filter.py gff_file -c --sorted --tandem-filtration

    fga_gene_filter.py gff_file -c -m --sorted --tandem--filtration

Such two calls will produce the same output, as after removing non-contiguous sequences of gene ids, computing median of median
coordinates for the rest of gene ids, and then sorting them, will not changing the relative order of different gene ids,
as while we have a contiguous sequence of same gene ids, that sequence with rewritten coordinates will be located (relatively to other
gene ids) at the same place in the overall order.

As ``-h/--help`` options states:

    usage: fga_gene_filter.py [-h] [--good-gene-ids-file GOOD_GENE_IDS_FILE]
                          [--bad-gene-ids-file BAD_GENE_IDS_FILE]
                          [--tandem-filtration] [--sorted] [-m] [-c]
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
      --tandem-filtration   substitutes every tandem duplication of gene ids, with
                            just one copy
      --sorted              sorts all gene ids, that survived filtration, on
                            respective fragments
      -m, --median          rewrites each gene ids coordinates with median of
                            median among all same gene ids coordinates
      -c, --continuous      filters out all gene ids, that, if being sorted by bp
                            coordinates, don't form a contiguous sequence

the script expect a full path to gff formatted file, or standard input.

For correct performance of ``--tandem--duplication`` option, if used any of the ``-c/--continuous``, ``-m/--median`` options, one shall
 supply ``--sorted`` option as well.

[4]:https://github.com/sergey-aganezov-jr/fga/blob/master/src/gene/fga_gene_filter.py

## [fga_gene_mgra_i.py][5]

This script is designed to perform transition from gff formatted data into mgra suitable format. As mgra software
expects genomes in grimm/infercars format, this script, taking **gene<->number mapping** file and gff formatted files,
and outputs genomes, represented as a set of fragments, where each fragment is represented as a sorted sequence of gene
 families. For now scripts outputs only grimm based format. More on grimm format can be found
 [here](http://grimm.ucsd.edu/GRIMM/).
 For the purpose of translating gff formatted data into gene family notation, suitable for mgra, script requires
 a **gene<-> number mapping** file, which wil contain at least all gene ids, that will be found in gff formatted files
 (can contain more, but not less, as KeyError would be raised and script will fail during execution).

Scripts derives genome names from the base name of supplied gff formatted files. The algorithm for name deriving is
 as follows:

1. Get string, which precedes the last ``.`` char (omit the extension)
2. Get string, which in previous result, precedes the first ``_`` char, if any.


#### Usage

    fga_gene_mgra_i.py gene_number_mapping gff_file_one gff_file_two ...

This call will output the contents of gff\_file\_one and gff\_file\_two which will be suitable fo mgra. The output will
be as follows:

    >gff_file_one
    # content of gff_file_one

    >gff_file_two
    # content of gff_file_two

Script doesn't perform any additional work on supplied data, besides sorting it based on respective base pair
coordinates.

As ``-h/--help`` option states:

    usage: fga_gene_mgra_i.py [-h]
                          gene_number_mapping_file gff_files [gff_files ...]

    positional arguments:
      gene_number_mapping_file
                            full path to gene number mapping file. Format can be
                            seen in docs section,in gene section.
      gff_files             full path to gff files for different genomes

    optional arguments:
      -h, --help            show this help message and exit

Script expects **gene<->number mapping** file as it first argument and, as stated above, all gene ids, that would be
found in further analysed gff formatted data must be in the **gene<->number mapping** file.
Script expects at least one gff formatted file, as it second and further positional arguments. Script doesn't read from
standard input, as it tries to determine genome names from supplied gff formatted file names.

[5]:https://github.com/sergey-aganezov-jr/fga/blob/master/src/gene/fga_gene_mgra_i.py

# Assembly

This section covers all scripts that can be found in ``assembly`` python package in root ``src`` package. These scripts
are designed to process MGRA log files, and input raw files in conjunction, to output a readable information, about the
assembled info, that was gained by using genome rearrangement analysis tools, to perform gluing among different
 fragments in highly fragmented genomes.

MGRA log like file is expected to contain specifically formatted data. There shall be rows of two types:

 1. ``A=ABCD=Abcdefg``: this row represents so called genome aliases definition. It states three value, first of which
 will be used to determine in which genome particular gluing operations happened. Second four letter value determines
 the relationship between gene id mapped to some value and genome, it belongs to. Third value stands for the full genome
 name that will be obtained by processing basename of supplied gff formatted file.
 2. ``(12h, oo)x(11t, oo):{A}``: this row represents information about gluing operation in particular genome. One letter
    alias in brackets, while digital value prefixing ``h`` or ``t`` letter represents an end of outer-most synteny
     blocks of some fragments, that were glued together.

In MGRA log like file script expects to find genome aliases definition for at least all one letter aliases, that could
be found in gluing operation.

Output format:

    genome_1
        fragment_1 (+) <==> fragment_2 (-)

All rows that don't have tabulation in beginning, contain only genome name. All rows that start with tabulation
character, contain at least one fragment name from last mentioned genome (obtained from first column of
 gff formatted files). The ``+`` or ``-`` sign determines if fragment has to be oriented positively, or
  is it has to be reversed (``-`` sign can be found only in assembled chains of length at
 least 2). If any number of fragments in one row is followed by ``<==>``, this means that previously mentioned
 sequence of fragments was glued together with a fragment, that follows the ``<==>``. In the example above
 ``fragment_1`` oriented naturally, was glued with ``fragment_2`` in reversed position. In other terms, if each
 fragment can be represented as a pair of ``start`` and ``end`` values, example above shows that ``end`` of
 ``fragment_1`` was glued together with ``end`` of ``fragment_2``.

## Usage

    fga_assmebly_fgr.py mgra_log_file gene_number_mapping_file gff_file_one gff_file_two ...

This call will analyse mgra log like file, obtain information about which fragments were glued with which fragments,
and in which genomes, than, using the information from gene_number_mapping file, it will translate mgra integer based
synteny blocks names, into the gene_ids based name. Using the knowledge if which gene ids were used in gluing operations,
script determines from which fragments in which genomes were those gene ids. Script double checks, that those gene ids
 indeed represents the outermost gene ids on respective fragments and that there are no cases of fragment being glued to
 another end of itself, or that some fragment took place in more than 2 gluing operations, as it has only two
 extremities, thus only two possible options to participate in gluing operations. After determining which fragments
 were assembled together script reports them as glued together chains.

As ``-h/--help`` option states:

    usage: fga_assembly_fgr.py [-h] [-c]
                           mgral_file gene_num_mapping_file gff_files
                           [gff_files ...]

    positional arguments:
      mgral_file            full path to mgra like log file with information about
                            glued vertices
      gene_num_mapping_file
                            full path to gene mapping file
      gff_files             full path to gff formatted info files

    optional arguments:
      -h, --help            show this help message and exit
      -c, --chains          permits output for assembled fragments only

Script expects exactly one **mgra log like** file, one **gene<->number mapping** file and at least one **gff formatted**
 file. Scripts outputs each supplied in gff format genome as a set of fragments, that are obtained from first column in
 those gff formatted files.

Just like the scripts from ``gene`` section, ``fga_assembly_fgr.py`` script, while reading genome information from gff
formatted files, determines the genome name from the gff formatted file name, using the same algorithm for determining
 it. Its important to notice, that this name must be the same, as the third value in the genome aliases definition row
 in mgra log like file.

Using the ``-c`` option one can limit the output of each genome to only those fragments, that were a part of at least
 one gluing.





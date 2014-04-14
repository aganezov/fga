# Documentation

This is a documentation section for **fga** project software, which will help you navigate around source code,
understand basic principals of what which script does, and hopefully answer all you questions.

## Purpose

The **FGA** project was undertaken to study evolutionary genomics in several mosquito genomes form South Africa.

As the data had to be prepared and analyzed, several software scripts for those purposes were created. All scripts
were written to be used as command line stand alone tools.

In general, all software is divided into following categories:

1. **orth** - this set of scripts deals with data, that contains information about orthology relationship among genomes
    (gene families). More on particular script usage and data format in [orth documentation][1]
2. **gene** - this set of scripts deals with data, that contains information about coding exons and fragments in each
    genome. More on particular script usage and data format in [gene documentation][2]
3. **assembly** - this set of scripts deals with data, that contains information about assembled statistics, among genomes
    More on particular script usage and data format in [assembly documentation][3]

Each set is contained in a respectively named python package. All scripts are written to be invoked as standalone command
line tools. All scripts write to standard output, thus can be used as unix pipes sources. Some scripts can read from standard
input only, if needed, thus can be used in any position in unix pipes.

All scripts are written with support of python3.3+. Each script contains ``#! /usr/bin/env python3`` directive, thus can be
invoked without using ``python3 script_name.py`` prefix.

As each script was meant to be standalone and not dependent on any other package, library (except standard), etc., code
contains some duplication in terms of support functions.

## Usage

Please find usage examples for each particular group and each individual scripts in respective documentation packages.

[1]:https://github.com/sergey-aganezov-jr/fga/tree/master/docs/orth
[2]:https://github.com/sergey-aganezov-jr/fga/tree/master/docs/gene
[3]:https://github.com/sergey-aganezov-jr/fga/tree/master/docs/assembly
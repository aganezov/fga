#fgr

this is a quick overview on how to use, and what to expect from fga_assembly_fgr.py script.

Script takes as input several data files with several genome information and outputs each genome as a set of fragments,
where some oos those are represented as oriented sequence of smaller fragments (result of assembly).

Using **-c** command line option when invoking script will limit its output to only those fragments in each genomes,
 that were gained, by gluing at least two fragments in original data.

As **-h** options shows, scripts awaits for at least three positional arguments:

1. mgral file -- MGRA-log-like file, which contains information about fusions, that are performed in different genomes
as well as genomes definitions

        \# genome definition part
        R=AARA=arabiensis

        \# gluing statistics
        (1604t,oo)x(360h,oo):{R}  (1691h,oo)x(4304t,oo):{R}  (3197h,oo)x(521h,oo):{R}  (2001t,oo)x(3889h,oo):{R}
        (4070h,oo)x(597h,oo):{R}  (173t,oo)x(6061h,oo):{R}  (4375h,oo)x(768h,oo):{R}  (1182h,oo)x(1767t,oo):{R}
        (127h,oo)x(963t,oo):{R}  (1950t,oo)x(982t,oo):{R}

        * R stands for one letter alias for arabiensis organism
        * AARA stands for gene_id prefix, that is used in gff file for arabiensis
        * arabiensis stands for verbose name for this particular creature

    for more examples look at gluing.txt file in current directory

2. gene_id -\> number mapping file - this file contains simple pairs of gene_ids and numbers, that those gene were
 mapped to for MGRA input purposes. This file doesn't allow any other string, but such mapping pairs.

        AARA000009 4145

    for more examples look at gene_num_map.txx file in current directory

3. gff files for each particular genome. More on gff format [here](http://www.ensembl.org/info/website/upload/gff.html)

Algorithm **checks**, that the number of genome definitions in mgral file equals to the number of genomes, that are
extracted from one-letter aliases in gluing statistics.

Algorithm **relies**, that in genome definition gene_id prefixes are of 4 letter length, as in gff files, and gene number
mapping file it parses each gene (to determine which genome it belong to) by it's 4 letter prefix.

Algorithm **checks**, that in gluing statistics each breakpoint graph vertex (represented by a particular number), does
 reflect a starting or ending point of genome fragment.

Algortihm **chekcs**, that in gluing statistics each breakpoint graph vertex (represented by a particular number), does
reflect a respective gene end, as determined by strand, gene is located on, and "t"|"h" attribute for particular vertex.







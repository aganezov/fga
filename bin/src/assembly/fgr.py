# -*- coding: utf-8 -*-
from collections import defaultdict
import os
import sys
import re

__author__ = 'Sergey Aganezov Jr.'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'

"""
This script shall be invoked as a command line tool.
This script is written with python3.4+ support.
This script accepts only command line arguments.

Idea:
    The idea is to take MGRA-log-kind-of file (MGRAL), several gff files, and gene mapping file,
    to report the fragments, that were glued together.
    Script analyses MGRAL file, identifying fusion operations and vertices they affect, after that script maps MGRA
    BG vertices representation with actual gene ids they correspond to. After that gff files are processed in order
    to map respective genes to fragments they are located on, perform several checks, and report which fragments were
    glued together.

    Example:
        MGRAL file states:
                R=AARA=Arabiensis
                ...
                (1t, oo)x(oo, 2t):{R}

        Gene mapping file states:
                AARA0078 1
                AARA0098 2
        gff file for Arabiensis states:
            K2345 ... AARA0078  .. 1 5 +
            K2356 ... AARA0098  .. 1 10 +

        script will report:
            in Arabiensis:
                K2345(-) <==> K2356(+)

    double checks:
        script checks that genes, gluing is performed on, are left/right most genes on fragments
        script checks that genes extremities, gluing is performed on, represent the end of the segment
        script check that genes, gluing is performed on, are on different fragments


Input:


Output:

"""

# yeap, evil
vrregex = re.compile("\((?:(?P<vertex1>\d+(?:h|t))|oo),(?:(?P<vertex2>\d+(?:h|t))|oo)\)x\((?:(?P<vertex3>\d+(?:h|t))|oo),(?:(?P<vertex4>\d+(?:h|t))|oo)\):\{(?P<genome_alias>\w)\}")
gnsregex = re.compile("(?P<ola>\w)=(?P<gkn>\w+)=(?P<gn>\w+)")


def retrieve_gene_number_mapping_from_file(gene_mapping_file):
    result = defaultdict(dict)
    with open(gene_mapping_file, "r") as source:
        for line in source:
            line = line.strip()
            gene_id, number = line.split()
            number = int(number)
            key = gene_id[:4]
            result[key][number] = gene_id
    return result


def remove_tandem_duplications(fragment_content):
    if len(fragment_content) == 1:
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
    assert set(x[0] for x in result_genome_aliases) == set(result_gluing_lists.keys())
    return result_genome_aliases, result_gluing_lists


def retrieve_genome_from_file(gff_file):
    result = defaultdict(list)
    with open(gff_file, "r") as source:
        for line in source:
            line = line.strip()
            data = line.split("\t")
            fragment = data[0]
            start, end = int(data[3]), int(data[4])
            strand = data[6]
            gene_id = data[10]
            result[fragment].append((gene_id, start, end, strand))
    for fragment in result:
        result[fragment] = sorted(result[fragment], key=lambda x: (x[1], x[2]))
        result[fragment] = remove_tandem_duplications(result[fragment])
    return result


def find_fragment_by_extremity(genome, gene_id, vertex):
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


def main(mgral_file, gene_mapping_file, gff_files):
    gene_number_mapping = retrieve_gene_number_mapping_from_file(gene_mapping_file=gene_mapping_file)
    genome_aliases, gluing_results = retrieve_glued_vertices_from_file(mgral_file=mgral_file)
    genomes = {}
    result = defaultdict(list)
    for gff_file in gff_files:
        genome_name = os.path.basename(gff_file).split(".")[0]
        genomes[genome_name] = retrieve_genome_from_file(gff_file)
    for genome_one_letter_alias in gluing_results:
        genome_name = genome_aliases[genome_one_letter_alias][1]
        genome_key = genome_aliases[genome_one_letter_alias][0]
        for v1, v2 in gluing_results[genome_one_letter_alias]:
            fragment1, d1 = find_fragment_by_extremity(genomes[genome_name], gene_number_mapping[genome_key][int(v1[:-1])], v1)
            fragment2, d2 = find_fragment_by_extremity(genomes[genome_name], gene_number_mapping[genome_key][int(v2[:-1])], v2)
            # assert fragment1 != fragment2
            result[genome_name].append((fragment1, fragment2, d1, d2))
    for genome in result:
        print(genome)
        gluing_statistics = result[genome]
        for f1, f2, d1, d2 in gluing_statistics:
            if f1 == f2:
                print("\t", f1, "!!! (circularized)")
            if d1 == 1:
                if d2 == -1:
                    print("\t", f1, "(+) <==>", f2, "(+)")
                else:
                    print("\t", f1, "(+) <==>", f2, "(-)")
            else:
                if d2 == -1:
                    print("\t", f1, "(-) <==>", f2, "(+)")
                else:
                    print("\t", f2, "(+) <==>", f1, "(-)")




if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    main(cmd_args[0], cmd_args[1], cmd_args[2:])
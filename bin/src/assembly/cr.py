# -*- coding: utf-8 -*-
from collections import defaultdict
import sys
import re

__author__ = "Sergey Aganezov"
__email__ = "aganezov@gwu.edu"
__status__ = "develop"

"""
"""

gnregex = re.compile("(?P<genome>\w+)")
gregex = re.compile("\t(?P<f1>\w+) \((?P<d1>\+|\-)\) \<\=\=\> (?P<f2>\w+) \((?P<d2>\+|\-)\)")

def main(input_file=None):
    if input_file is not None:
        source = open(input_file, "r")
    else:
        source = sys.stdin

    data = source.readlines()
    data = list(map(lambda line: line.rstrip(), data))
    parsed_data = defaultdict(list)
    current_genome = None
    for line in data:
        if gnregex.match(line):
            current_genome = gnregex.search(line).groupdict()["genome"]
        elif gregex.match(line):
            search_result = gregex.search(line)
            fragment_1, fragment_2 = search_result.groupdict()["f1"], search_result.groupdict()["f2"]
            d1, d2 = search_result.groupdict()["d2"], search_result.groupdict()["d1"]
    if not source is sys.stdin:
        source.close()


if __name__ == "__main__":
    cmd_args = sys.argv[:1]

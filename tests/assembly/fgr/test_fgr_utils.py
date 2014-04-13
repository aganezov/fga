# -*- coding: utf-8 -*-
from src.assembly import fga_assembly_fgr

__author__ = 'Sergey Aganezov Jr.'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'

import unittest


class fgrTandemDuplicationsTestCase(unittest.TestCase):
    def setUp(self):
        self.genome = {
            "fragment_1": [("gene_1", 1, 5, "+"), ("gene_1", 6, 8, "+"), ("gene_2", 9, 12, "-")],
            "fragment_2": [("gene_2", 5, 10, "-"), ("gene_1", 11, 20, "+"), ("gene_1", 23, 25, "+")],
            "fragment_3": [("gene_4", 1, 10, "-"), ("gene_3", 13, 15, "-"), ("gene_2", 17, 20, "-")],
            "fragment_4": [("gene_1", 1, 2, "-"), ("gene_2", 3, 4, "-"), ("gene_3", 5, 6, "-")],
            "fragment_5": [("gene_1", 1, 2, "-"), ("gene_2", 3, 4, "-"), ("gene_2", 5, 6, "+"), ("gene_3", 7, 8, "+")]
        }

    def test_collapse_tandem_duplication_at_fragment_start(self):
        duplication_at_start = fga_assembly_fgr.colapse_tandem_duplications(self.genome["fragment_1"])
        self.assertEqual(len(duplication_at_start), 2)
        self.assertListEqual(duplication_at_start, self.genome["fragment_1"][::2])

    def test_collapse_tandem_duplication_at_fragment_end(self):
        duplication_at_end = fga_assembly_fgr.colapse_tandem_duplications(self.genome["fragment_2"])
        self.assertEqual(len(duplication_at_end), 2)
        self.assertListEqual(duplication_at_end, self.genome["fragment_2"][:-1])

    def test_collapse_tandem_duplication_at_fragment_center(self):
        duplication_at_center = fga_assembly_fgr.colapse_tandem_duplications(self.genome["fragment_5"])
        self.assertEqual(len(duplication_at_center), 3)
        self.assertListEqual(duplication_at_center, self.genome["fragment_5"][:2] + self.genome["fragment_5"][-1:])

    def test_collapse_tandem_duplication_no_occurrence(self):
        no_occurrence = fga_assembly_fgr.colapse_tandem_duplications(self.genome["fragment_3"])
        self.assertListEqual(no_occurrence, self.genome["fragment_3"])
        no_occurrence = fga_assembly_fgr.colapse_tandem_duplications(self.genome["fragment_4"])
        self.assertListEqual(no_occurrence, self.genome["fragment_4"])



if __name__ == '__main__':
    unittest.main()

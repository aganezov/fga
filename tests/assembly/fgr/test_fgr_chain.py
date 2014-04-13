# -*- coding: utf-8 -*-
from collections import defaultdict

__author__ = 'Sergey Aganezov'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'
import unittest
from src.assembly import fga_assembly_fgr


class fgrConnectionCreationTestCase(unittest.TestCase):
    def setUp(self):
        self.storage = defaultdict(lambda: defaultdict(list))

    def test_forward_forward_case(self):
        fragment_1, fragment_2 = "f1", "f2"
        d1, d2 = 1, -1
        fga_assembly_fgr.connection_creation(f1=fragment_1, f2=fragment_2, d1=d1, d2=d2, storage=self.storage)
        self.assertEqual(len(self.storage), 2)
        self.assertTrue(fragment_1 in self.storage)
        self.assertTrue(fragment_2 in self.storage)
        self.assertTrue("h" in self.storage[fragment_1])
        self.assertTrue("t" in self.storage[fragment_2])
        self.assertListEqual(self.storage[fragment_1]["h"], [(fragment_2, "t")])
        self.assertListEqual(self.storage[fragment_2]["t"], [(fragment_1, "h")])

    def test_forward_backwards_case(self):
        fragment_1, fragment_2 = "f1", "f2"
        d1, d2 = 1, 1
        fga_assembly_fgr.connection_creation(f1=fragment_1, f2=fragment_2, d1=d1, d2=d2, storage=self.storage)
        self.assertEqual(len(self.storage), 2)
        self.assertTrue(fragment_1 in self.storage)
        self.assertTrue(fragment_2 in self.storage)
        self.assertTrue("h" in self.storage[fragment_1])
        self.assertTrue("h" in self.storage[fragment_2])
        self.assertListEqual(self.storage[fragment_1]["h"], [(fragment_2, "h")])
        self.assertListEqual(self.storage[fragment_2]["h"], [(fragment_1, "h")])

    def test_backwards_forward(self):
        fragment_1, fragment_2 = "f1", "f2"
        d1, d2 = -1, 1
        fga_assembly_fgr.connection_creation(f1=fragment_1, f2=fragment_2, d1=d1, d2=d2, storage=self.storage)
        self.assertEqual(len(self.storage), 2)
        self.assertTrue(fragment_1 in self.storage)
        self.assertTrue(fragment_2 in self.storage)
        self.assertTrue("t" in self.storage[fragment_1])
        self.assertTrue("h" in self.storage[fragment_2])
        self.assertListEqual(self.storage[fragment_1]["t"], [(fragment_2, "h")])
        self.assertListEqual(self.storage[fragment_2]["h"], [(fragment_1, "t")])

    def test_backwards_backwards_case(self):
        fragment_1, fragment_2 = "f1", "f2"
        d1, d2 = -1, -1
        fga_assembly_fgr.connection_creation(f1=fragment_1, f2=fragment_2, d1=d1, d2=d2, storage=self.storage)
        self.assertEqual(len(self.storage), 2)
        self.assertTrue(fragment_1 in self.storage)
        self.assertTrue(fragment_2 in self.storage)
        self.assertTrue("t" in self.storage[fragment_1])
        self.assertTrue("t" in self.storage[fragment_2])
        self.assertListEqual(self.storage[fragment_1]["t"], [(fragment_2, "t")])
        self.assertListEqual(self.storage[fragment_2]["t"], [(fragment_1, "t")])


class fgrConnectionConstructionTestCase(unittest.TestCase):
    def setUp(self):
        self.genome = {
            "fragment_1": [("gene_1", 1, 5, "+"), ("gene_2", 6, 8, "+"), ("gene_3", 9, 12, "-")],
            "fragment_2": [("gene_4", 5, 10, "-"), ("gene_5", 11, 20, "+"), ("gene_6", 23, 25, "+")],
            "fragment_3": [("gene_7", 1, 10, "-"), ("gene_8", 13, 15, "-"), ("gene_9", 17, 20, "-")],
            "fragment_4": [("gene_10", 1, 2, "-"), ("gene_11", 3, 4, "-"), ("gene_12", 5, 6, "-")],
            "fragment_5": [("gene_13", 1, 2, "-"), ("gene_14", 3, 4, "-"), ("gene_15", 5, 6, "+"),
                           ("gene_16", 7, 8, "+")],
            "fragment_6": [("gene_17", 1, 2, "-"), ("gene_18", 3, 4, "-"), ("gene_19", 5, 6, "-")],
        }
        self.connection_info_storage = {"fragment_1": {"t": [("fragment_2", "t")],
                                                       "h": [("fragment_3", "t")]},
                                        "fragment_2": {"t": [("fragment_1", "t")]},
                                        "fragment_3": {"t": [("fragment_1", "h")]},
                                        "fragment_5": {"h": [("fragment_6", "t")]},
                                        "fragment_6": {"t": [("fragment_5", "h")]}
        }
        self.connection_info_storage = defaultdict(lambda: defaultdict(list))
        self.connection_info_storage["fragment_1"]["t"].append(("fragment_2", "t"))
        self.connection_info_storage["fragment_1"]["h"].append(("fragment_3", "t"))
        self.connection_info_storage["fragment_2"]["t"].append(("fragment_1", "t"))
        self.connection_info_storage["fragment_3"]["t"].append(("fragment_1", "h"))
        self.connection_info_storage["fragment_5"]["h"].append(("fragment_6", "t"))
        self.connection_info_storage["fragment_6"]["t"].append(("fragment_5", "h"))

        self.visited = {fragment: False for fragment in self.genome}

    def test_chain_construction_no_connections(self):
        result_chain = fga_assembly_fgr.chain_construction("fragment_4", None, storage=self.connection_info_storage,
                                                           visited=self.visited)
        self.assertEqual(len(result_chain), 1)
        self.assertListEqual([("fragment_4", "+")], result_chain)

    def test_chain_construction_single_gluing(self):
        result_chain = fga_assembly_fgr.chain_construction("fragment_5", "h", storage=self.connection_info_storage,
                                                           visited=self.visited)
        self.assertEqual(len(result_chain), 2)
        self.assertListEqual([("fragment_5", "+"), ("fragment_6", "+")], result_chain)

        self.visited = {fragment: False for fragment in self.genome}
        result_chain = fga_assembly_fgr.chain_construction("fragment_6", "t", storage=self.connection_info_storage,
                                                           visited=self.visited)
        self.assertEqual(len(result_chain), 2)
        self.assertListEqual([("fragment_6", "-"), ("fragment_5", "-")], result_chain)

    def test_chain_construction_fails_on_non_starting_fragment(self):
        with self.assertRaises(AssertionError):
            fga_assembly_fgr.chain_construction("fragment_1", "h", storage=self.connection_info_storage,
                                                visited=self.visited)
        with self.assertRaises(AssertionError):
            fga_assembly_fgr.chain_construction("fragment_1", None, storage=self.connection_info_storage,
                                                visited=self.visited)

    def test_chain_construction_multiple_gluings(self):
        result_chain = fga_assembly_fgr.chain_construction("fragment_2", "t", storage=self.connection_info_storage,
                                                           visited=self.visited)
        self.assertEqual(len(result_chain), 3)
        self.assertListEqual([("fragment_2", "-"), ("fragment_1", "+"), ("fragment_3", "+")], result_chain)

        self.visited = {fragment: False for fragment in self.genome}
        result_chain = fga_assembly_fgr.chain_construction("fragment_3", "t", storage=self.connection_info_storage,
                                                           visited=self.visited)

        self.assertEqual(len(result_chain), 3)
        self.assertListEqual([("fragment_3", "-"), ("fragment_1", "-"), ("fragment_2", "+")], result_chain)

    def test_visited_array_fill(self):
        def retrieve_only_true(dict):
            return {key: value for key, value in dict.items() if value}

        def retrieve_only_false(dict):
            return {key: value for key, value in dict.items() if not value}

        self.assertEqual(len(retrieve_only_false(self.visited)), len(self.genome))
        self.assertEqual(len(retrieve_only_true(self.visited)), 0)

        fga_assembly_fgr.chain_construction("fragment_4", None, storage=self.connection_info_storage,
                                            visited=self.visited)

        self.assertTrue(self.visited["fragment_4"])
        self.assertTrue(len(retrieve_only_false(self.visited)), len(self.genome) - 1)

        result_chain = fga_assembly_fgr.chain_construction("fragment_2", "t", storage=self.connection_info_storage,
                                                           visited=self.visited)
        for fragment, d in result_chain:
            self.assertTrue(self.visited[fragment])



if __name__ == '__main__':
    unittest.main()

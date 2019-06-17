import unittest

from phylogenetics import UltrametricTree


class test_UltrametricTree(unittest.TestCase):

    a = UltrametricTree(1)
    b = UltrametricTree(2)
    c = UltrametricTree(3)
    d = UltrametricTree(4)
    e = UltrametricTree(-1, 2.0, b, c)
    f = UltrametricTree(-1, 1.5, d, a)
    final_tree = UltrametricTree(-1, 4.5, e, f)

    def test_get_cluster(self):
        self.assertEqual(self.e.get_cluster(), [2, 3])
        self.assertEqual(self.final_tree.get_cluster(), [2, 3, 4, 1])

        print('>> Passed test_get_cluster()')

    def test_size(self):
        self.assertEqual(self.e.size(), (1, 2))
        self.assertEqual(self.final_tree.size(), (3, 4))

        print('>> Passed test_size()')

    def test_exists_leaf(self):
        self.assertFalse(self.e.exists_leaf(1))
        self.assertTrue(self.final_tree.exists_leaf(3))

        print('>> Passed test_exists_leaf()')

    def test_common_ancestor(self):
        self.assertEqual(self.final_tree.common_ancestor(2, 3), self.e)
        self.assertEqual(self.final_tree.common_ancestor(1, 3), self.final_tree)

        print('>> Passed test_common_ancestor()')

    def test_distance_leaves(self):
        self.assertEqual(self.final_tree.distance_leaves(2, 3), 4.0)
        self.assertEqual(self.final_tree.distance_leaves(1, 3), 9.0)

        print('>> Passed test_distance_leaves()')

if __name__ == '__main__':
    unittest.main()

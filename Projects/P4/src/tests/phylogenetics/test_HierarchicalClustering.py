import unittest

from phylogenetics import HierarchicalClustering, NumMatrix


class test_HierarchicalClustering(unittest.TestCase):

    m = NumMatrix(5, 5)
    m.set_value(0, 1, 2)
    m.set_value(0, 2, 5)
    m.set_value(0, 3, 7)
    m.set_value(0, 4, 9)
    m.set_value(1, 2, 4)
    m.set_value(1, 3, 6)
    m.set_value(1, 4, 7)
    m.set_value(2, 3, 4)
    m.set_value(2, 4, 6)
    m.set_value(3, 4, 3)
    hc = HierarchicalClustering(m)

    def test_execute_clustering(self):
        tree = self.hc.execute_clustering()
        self.assertTrue(tree.exists_leaf(1))
        self.assertTrue(tree.exists_leaf(4))
        self.assertTrue(tree.exists_leaf(3))
        self.assertEqual(tree.distance_leaves(1, 3), 6.5)

        print('>> Passed test_execute_clustering()')


if __name__ == '__main__':
    unittest.main()

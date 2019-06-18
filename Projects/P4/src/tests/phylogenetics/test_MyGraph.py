import unittest

from ...phylogenetics import MyGraph, NumMatrix


class test_MyGraph(unittest.TestCase):

    gr = MyGraph()
    gr.add_node(1)
    gr.add_node(2)
    gr.add_node(3)
    gr.add_node(4)
    gr.add_edge(1, 2)
    gr.add_edge(2, 3)
    gr.add_edge(3, 2)
    gr.add_edge(3, 4)
    gr.add_edge(4, 2)

    def test_basics(self):
        self.assertEqual([1, 2, 3, 4], self.gr.get_nodes())
        self.assertEqual(
            [(1, 2), (2, 3), (3, 2), (3, 4), (4, 2)],
            self.gr.get_edges())

        self.assertEqual(self.gr.size(), (4, 5))

        print('>> Passed MyGraph::test_basics()')

    def test_neighborhood(self):
        self.assertEqual(self.gr.get_successors(2), [3])
        self.assertEqual(self.gr.get_predecessors(2), [1, 3, 4])
        self.assertEqual(self.gr.get_adjacents(2), [1, 3, 4])

        print('>> Passed test_neighborhood()')

    def test_degree(self):
        self.assertEqual(self.gr.out_degree(2), 1)
        self.assertEqual(self.gr.in_degree(2), 3)
        self.assertEqual(self.gr.degree(2), 4)
        self.assertEqual(self.gr.all_degrees()[2], 4)
        self.assertEqual(self.gr.highest_degrees()[0], 2)

        print('>> Passed test_degree()')

    def test_topological(self):
        self.assertEqual(self.gr.mean_degree(), 2.5)
        self.assertEqual(self.gr.prob_degree(),
                         {1: 0.25, 4: 0.25, 3: 0.25, 2: 0.25})

        print('>> Passed test_topological()')

    def test_searches(self):
        gr2 = MyGraph({1: [2, 3, 4], 2: [5, 6], 3: [6, 8],
                       4: [8], 5: [7], 6: [], 7: [], 8: []})

        self.assertEqual(gr2.reachable_bfs(1), [2, 3, 4, 5, 6, 8, 7])
        self.assertEqual(gr2.reachable_dfs(1), [2, 5, 7, 6, 3, 8, 4])
        self.assertEqual(gr2.distance(1, 7), 3)
        self.assertEqual(gr2.shortest_path(1, 7), [1, 2, 5, 7])
        self.assertEqual(gr2.distance(1, 8), 2)
        self.assertEqual(gr2.shortest_path(1, 8), [1, 3, 8])
        self.assertEqual(gr2.distance(6, 1), None)
        self.assertEqual(gr2.shortest_path(6, 1), None)
        self.assertEqual(gr2.mean_distances(), 2.5)

        print('>> Passed test_searches()')

    def test_clustering(self):
        self.assertEqual(self.gr.clustering_coef(1), 0)
        self.assertEqual(self.gr.clustering_coef(3), 1)
        self.assertEqual(self.gr.all_clustering_coefs(),
                         {1: 0.0, 2: 0.3333333333333333, 3: 1.0, 4: 1.0})
        self.assertEqual(self.gr.mean_clustering_coef(), 0.5833333333333333)

        print('>> Passed test_clustering()')

    def test_create_from_num_matrix(self):
        mat = NumMatrix(4, 4)
        mat.set_value(0, 1, 1)
        mat.set_value(1, 2, 5)
        mat.set_value(1, 3, 7)

        gr3 = MyGraph.create_from_num_matrix(mat, 6)
        self.assertEqual(
            gr3.get_edges(),
            [(1, 0), (2, 0), (2, 1), (3, 0), (3, 2)])
        self.assertTrue((1, 3) not in gr3.get_edges())

        print('>> Passed test_create_from_num_matrix()')


if __name__ == '__main__':
    unittest.main()

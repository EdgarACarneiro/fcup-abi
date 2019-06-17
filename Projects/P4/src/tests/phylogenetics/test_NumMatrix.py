import unittest

from phylogenetics import NumMatrix


class test_NumMatrix(unittest.TestCase):

    mat = NumMatrix(4, 4)

    def test_basics(self):
        self.assertEqual(self.mat[1], [0.0, 0.0, 0.0, 0.0])
        self.assertEqual(self.mat.num_rows(), 4)
        self.assertEqual(self.mat.num_cols(), 4)

        self.mat.set_value(2, 3, 7.0)
        self.assertEqual(self.mat.get_value(2, 3), 7.0)
        self.assertEqual(self.mat.get_value(3, 2), 7.0)

        self.mat.add_row([0, 0, 0, 0])
        self.assertEqual(self.mat.num_rows(), 5)

        self.mat.add_col([1, 1, 1, 1, 1])
        self.assertEqual(self.mat.num_cols(), 5)
        self.assertEqual(self.mat.get_value(4, 4), 1)

        self.mat.remove_row(4)
        self.assertEqual(self.mat.num_rows(), 4)

        self.mat.remove_col(4)
        self.assertEqual(self.mat.num_cols(), 4)

        print('>> Passed NumMatrix::test_basics()')

    def test_min_dist_indexes(self):
        self.mat.set_value(3, 2, -3.0)
        self.assertEqual(self.mat.min_dist_indexes(), (3, 2))

        print('>> Passed min_dist_indexes()')

    def test_copy(self):
        copy_mat = self.mat.copy()

        self.assertFalse(self.mat == copy_mat)
        self.assertTrue(self.mat.mat == copy_mat.mat)

        print('>> Passed test_copy()')

if __name__ == '__main__':
    unittest.main()

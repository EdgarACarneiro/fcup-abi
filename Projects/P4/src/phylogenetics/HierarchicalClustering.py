from . import UltrametricTree, NumMatrix


class HierarchicalClustering:
    """Implements a general purpose agglomerative hierarchical clustering"""

    dists_mat: NumMatrix

    def __init__(self, dists_mat):
        self.dists_mat = dists_mat

    def execute_clustering(self):
        """Runs the general purpose agglomerative hierarchical clustering algorithm
        and returns a UltrametricTree"""
        # initialize the trees
        trees = []
        for i in range(self.dists_mat.num_rows()):
            trees.append(UltrametricTree(i))  # All leafs

        # make a copy of the distance matrix to change it
        tableDist = self.dists_mat.copy()

        # iterations
        for k in range(self.dists_mat.num_rows(), 1, -1):
            mins = tableDist.min_dist_indexes()
            i, j = mins[0], mins[1]

            # create a new tree joining the clusters
            n = UltrametricTree(
                -1, tableDist.get_value(i, j) / 2.0, trees[i], trees[j])

            if k > 2:
                # remove trees being joined from the list
                ti = trees.pop(i)
                tj = trees.pop(j)
                dists = []

                # calculate the distance for the new cluster
                for x in range(tableDist.num_rows()):
                    if x != i and x != j:
                        si = len(ti.get_cluster())
                        sj = len(tj.get_cluster())

                        # use the weighted average to calculate the distances
                        # between the clusters
                        dists.append(
                            (si * tableDist.get_value(i, x) +
                                sj * tableDist.get_value(j, x)) / (si+sj))

                # update the matrix:
                tableDist.remove_row(i)
                tableDist.remove_col(i)
                tableDist.remove_row(j)
                tableDist.remove_col(j)

                tableDist.add_row(dists)
                tableDist.add_col([0] * (len(dists)+1))
                trees.append(n)

            else:
                return n

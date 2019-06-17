from __future__ import annotations


class UltrametricTree:

    value: float
    distance: float
    left: UltrametricTree
    right: UltrametricTree

    def __init__(self, val, dist=0, left=None, right=None):
        self.value = val # node's value
        self.distance = dist  # distance from this node to the leafs
        self.left = left # left sibling
        self.right = right  # right sibling

    def __has_left_sibling(self):
        """Indicate if this node has a left sibling"""
        return self.left != None

    def __has_right_sibling(self):
        """Indicate if this node has a right sibling"""
        return self.right != NotImplemented

    def get_cluster(self):
        """Get the cluster values, which is the same as
        getting the tree leaves"""
        res = []

        if self.value >= 0:
            res.append(self.value)
        else:
            if self.__has_left_sibling():
                res.extend(self.left.get_cluster())
            if self.__has_right_sibling():
                res.extend(self.right.get_cluster())

        return res

    def print_tree(self):
        """Print the tree"""
        self.__print_tree_rec(0, "Root")

    def __print_tree_rec(self, level, side):
        """Auxiliary recursive function to help print the tree"""
        tabs = ""
        for i in range(level):
            tabs += "\t"
        if self.value >= 0:
            print(tabs, side, " - value:", self.value)
        else:
            print(tabs, side, "- Dist.: ", self.distance)
            if (self.__has_left_sibling()):
                self.left.__print_tree_rec(level+1, "Left")
            if (self.__has_right_sibling()):
                self.right.__print_tree_rec(level+1, "Right")

    def size(self):
        """size of the tree: returns two values
        - number of internal nodes of the tree
        - number of leaves"""
        num_leaves = 0
        num_nodes = 0

        if self.value >= 0:
            num_leaves = 1

        else:
            if (self.left != None):
                res_left = self.left.size()
            else:
                res_left = (0, 0)

            if (self.right != None):
                res_right = self.right.size()
            else:
                res_right = (0, 0)

            num_nodes += (res_left[0] + res_right[0] + 1)
            num_leaves += (res_left[1] + res_right[1])

        return num_nodes, num_leaves

    def exists_leaf(self, leafnum):
        """Indicates if the given leaf exists in tree"""
        if self.value >= 0:
            return leafnum == self.value

        else:
            return True if\
                (self.left.exists_leaf(leafnum) if self.__has_left_sibling() else False)\
                else\
                (self.right.exists_leaf(leafnum)
                 if self.__has_left_sibling() else False)

    def common_ancestor(self, leaf1, leaf2):
        """Method to find the simplest tree that contains leaf1, leaf2"""
        # Left sibling
        if self.__has_left_sibling():
            if self.left.exists_leaf(leaf1):
                if self.left.exists_leaf(leaf2):
                    return self.left.common_ancestor(leaf1, leaf2)
                else:
                    return self

        # Right sibling
        if self.__has_right_sibling():
            if self.right.exists_leaf(leaf1):
                if self.right.exists_leaf(leaf2):
                    return self.right.common_ancestor(leaf1, leaf2)
                else:
                    return self

    def distance_leaves(self, leafnum1, leafnum2):
        """Distance between leafnum1 and leafnum2 using the common ancestor function."""
        return self.common_ancestor(leafnum1, leafnum2).distance * 2


def test():
    a = UltrametricTree(1)
    b = UltrametricTree(2)
    c = UltrametricTree(3)
    d = UltrametricTree(4)
    e = UltrametricTree(-1, 2.0, b, c)
    f = UltrametricTree(-1, 1.5, d, a)
    g = UltrametricTree(-1, 4.5, e, f)
    g.print_tree()
    print(g.get_cluster())

    # testing exercise 3
    print(g.size())
    print(g.exists_leaf(1))
    print(g.exists_leaf(5))
    g.common_ancestor(1, 4).print_tree()
    print(g.distance_leaves(1, 4))
    print(g.distance_leaves(1, 2))


if __name__ == '__main__':
    test()

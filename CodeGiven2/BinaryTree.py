class BinaryTree:

    def __init__(self, val, dist = 0, left = None, right = None):
        self.value = val
        self.distance = dist
        self.left = left
        self.right = right

    def has_left_sibling(self):
        return self.left != None

    def has_right_sibling(self):
        return self.right != NotImplemented

    def get_cluster(self):
        """Get the cluster values => get the tree leaves"""
        res = []

        if self.value >= 0:
            res.append(self.value)
        else:
            if self.has_left_sibling():
                res.extend(self.left.get_cluster())
            if self.has_right_sibling():
                res.extend(self.right.get_cluster())

        return res

    def print_tree(self):
        self.print_tree_rec(0, "Root")
    
    def print_tree_rec (self, level, side):
        tabs = ""
        for i in range(level): tabs += "\t"
        if self.value >= 0:
            print(tabs, side, " - value:", self.value)
        else:
            print(tabs, side, "- Dist.: ", self.distance)
            if (self.has_left_sibling()): 
                self.left.print_tree_rec(level+1, "Left")
            if (self.has_right_sibling()): 
                self.right.print_tree_rec(level+1, "Right")
     
    
    def size(self):
        ''' size of the tree: returns two values
        - number of internal nodes of the tree
        - number of leaves'''
        numleaves = 0
        numnodes = 0
        if self.value >= 0:
            numleaves = 1
        else: 
            if (self.left != None):
                resl = self.left.size()
            else: resl = (0,0)
            if (self.right != None):  
                resr = self.right.size() 
            else: resr = (0,0)
            numnodes += (resl[0] + resr[0] + 1)
            numleaves += (resl[1] + resr[1])
        return numnodes, numleaves

    def exists_leaf(self, leafnum) -> bool:
        """Finds leaf in tree"""

        if self.value >= 0:
            return leafnum == self.value
        else:
            return (self.left.exists_leaf(leafnum) if self.has_left_sibling() else False) or\
                    (self.right.exists_leaf(leafnum) if self.has_left_sibling() else False)
    
    def common_ancestor(self, leaf1, leaf2):
       ''' Return simplest tree that contains leaf1, leaf2'''
       pass
       # ...

    def distance_leaves(self, leafnum1, leafnum2):
       ''' distance between leafnum1 and leafnum2 using the common ancestor function.'''
       pass
       # ...

def test():              
    a = BinaryTree(1)
    b = BinaryTree(2)
    c = BinaryTree(3)
    d = BinaryTree(4)
    e = BinaryTree(-1, 2.0, b, c)
    f = BinaryTree(-1, 1.5, d, a)
    g = BinaryTree(-1, 4.5, e, f)
    g.print_tree()
    print(g.get_cluster())

    # testing exercise 3
    print(g.size())
    print(g.exists_leaf(1))
    print(g.exists_leaf(5))
    g.common_ancestor(1,4).print_tree()
    print(g.distance_leaves(1,4))
    print(g.distance_leaves(1,2))

if __name__ == '__main__':
    test()

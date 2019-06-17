class MyGraph:
    """Directed Graph represented as adjacency list using a dictionary.
    Keys are vertices.
    Values of the dictionary represent the list of adjacent vertices of the key node"""

    def __init__(self, g={}):
        """Constructor - takes dictionary to fill the graph as input; default is empty dictionary"""
        self.graph = g

    def print_graph(self):
        """Prints the content of the graph as adjacency list"""
        for v in self.graph.keys():
            print(v, " -> ", self.graph[v])

    # get basic info
    def get_nodes(self):
        return self.graph.keys()

    def get_edges(self):
        """Returns edges in the graph as a list of tuples (origin, destination)"""
        res = []
        for n in self.graph:
            res.extend([(n, n2) for n2 in self.graph[n]])

        return res

    def size(self):
        """Returns size of the graph : number of nodes, number of edges"""
        return len(self.graph), len(self.get_edges())

    # add nodes and edges
    def add_node(self, v):
        """Add a node to the graph; tests if node exists not adding if it does"""
        if (not v in self.graph.keys()):
            self.graph[v] = []

    def add_edge(self, o, d):
        """Add edge to the graph; if vertices do not exist, they are added to the graph"""
        self.add_node(o)
        self.add_node(d)
        if d not in self.graph[o]:
            self.graph[o].append(d)

    # successors, predecessors, adjacent nodes

    def get_successors(self, v):
        """Get the successors of the given node"""
        # needed to avoid list being overwritten of result of the function is used
        return list(self.graph[v])

    def get_predecessors(self, v):
        """Get the predecessors of the given node"""
        res = []
        for k in self.graph.keys():
            if v in self.graph[k]:
                res.append(k)
        return res

    def get_adjacents(self, v):
        """Get the adjacent nodes to the given node"""
        suc = self.get_successors(v)
        pred = self.get_predecessors(v)
        res = pred
        for p in suc:
            if p not in res:
                res.append(p)
        return res

    # degrees

    def out_degree(self, v):
        """Number of links the given node has leaving it"""
        return sum([1 for edge in self.get_edges() if edge[0] == v])

    def in_degree(self, v):
        """Number of links the given node has entering it"""
        return sum([1 for edge in self.get_edges() if edge[1] == v])

    def degree(self, v):
        """Number of links that are connected to the given node"""
        return sum([1 for edge in self.get_edges() if v in edge])

    def all_degrees(self, deg_type="inout"):
        """Computes the degree (of a given type) for all nodes.
        deg_type can be "in", "out", or "inout""""
        assert deg_type in ["inout", "in", "out"], "Invalid degree type"

        switcher = {
            "in": lambda v: self.in_degree(v),
            "out": lambda v: self.out_degree(v),
            "inout": lambda v: self.in_degree(v) + self.out_degree(v)
        }

        return {v: switcher[deg_type](v) for v in self.graph.keys()}

    def highest_degrees(self, all_deg=None, deg_type="inout", top=10):
        """Get the top 'top' nodes with the highest degree"""c
        if all_deg is None:
            all_deg = self.all_degrees(deg_type)

        ord_deg = sorted(list(all_deg.items()),
                         key=lambda x: x[1], reverse=True)
        return list(map(lambda x: x[0], ord_deg[:top]))

    # topological metrics over degrees

    def mean_degree(self, deg_type="inout"):
        """Average of all nodes degrees"""
        return sum(self.all_degrees(deg_type).values()) / self.size()[0]

    def prob_degree(self, deg_type="inout"):
        """Count the number of occurrences of each degree in the network
        and derive its frequencies"""
        probs = {}
        n_nodes = self.size()[0]

        for deg in self.all_degrees(deg_type).values():
            if (deg in probs):
                probs[deg] += 1 / n_nodes
            else:
                probs[deg] = 1 / n_nodes

        return probs

    # BFS and DFS searches

    def reachable_bfs(self, v):
        """Organize nodes in a breadth-first search style"""
        l = [v]   # list of nodes to be handled
        res = []  # list of nodes to return the result

        while len(l) > 0:
            node = l.pop(0)  # implements a queue: LILO

            if node != v:
                res.append(node)

            for elem in self.graph[node]:
                if elem not in res and elem not in l and elem != node:
                    l.append(elem)

        return res

    def reachable_dfs(self, v):
        """Organize nodes in a depth-first search style"""
        l = [v]
        res = []

        while len(l) > 0:
            node = l.pop(0)  # implements a stack:

            if node != v:
                res.append(node)

            s = 0
            for elem in self.graph[node]:
                if elem not in res and elem not in l:
                    l.insert(s, elem)
                    s += 1

        return res

    def distance(self, s, d):
        """Number of edges in teh shortest path between s and d"""
        if s == d:
            return 0

        l = [(s, 0)]
        visited = [s]

        while len(l) > 0:
            node, dist = l.pop(0)

            for elem in self.graph[node]:
                if elem == d:
                    return dist + 1
                elif elem not in visited:
                    l.append((elem, dist+1))
                    visited.append(elem)

        return None

    def shortest_path(self, s, d):
        """Path connecting the two nodes with the shortest length"""
        if s == d:
            return 0

        l = [(s, [])]
        visited = [s]

        while len(l) > 0:
            node, preds = l.pop(0)´

            for elem in self.graph[node]:
                if elem == d:
                    return preds+[node, elem]
                elif elem not in visited:
                    l.append((elem, preds+[node]))
                    visited.append(elem)

        return None

    def mean_distances(self):
        """Compute the average distance between the nodes in the network"""
        num_nodes = self.size()[0]

        return sum([self.distance(i, j)
                    for j in self.get_nodes()
                    for i in self.get_nodes()
                    if j > i]) / num_nodes

    # clustering

    def clustering_coef(self, v):
        """Indicates how connected are the given nodes on the network"""
        adjs = self.get_adjacents(v)
        if len(adjs) <= 1:
            return 0.0

        # calculate the number of links of the adjacent nodes
        # compare pairwisely if nodes in this list are connected between them
        ligs = sum([1 for j in adjs
                    for i in adjs
                    if i != j and i in self.get_adjacents(j)])

        return float(ligs)/(len(adjs)*(len(adjs)-1))

    def all_clustering_coefs(self):
        """Dictionary with all the clustering coefficient for all the nodes in the graph"""
        return {v: self.clustering_coef(v) for v in self.get_nodes()}

    def mean_clustering_coef(self):
        """Mean value of all clustering coefficient for all the nodes in the graph"""
        return self.all_clustering_coefs() / self.size()[0]


if __name__ == "__main__":
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
    gr.print_graph()
    print(gr.size())
    print()
    print(gr.get_successors(2))
    print(gr.get_predecessors(2))
    print(gr.get_adjacents(2))
    print()
    print(gr.in_degree(2))
    print(gr.out_degree(2))
    print(gr.degree(2))
    print()
    print(gr.all_degrees("inout"))
    print(gr.all_degrees("in"))
    print(gr.all_degrees("out"))
    print()
    #
    gr2 = MyGraph({1: [2, 3, 4], 2: [5, 6], 3: [6, 8],
                   4: [8], 5: [7], 6: [], 7: [], 8: []})
    print(gr2.reachable_bfs(1))
    print(gr2.reachable_dfs(1))
    print(gr2.distance(1, 7))
    print(gr2.shortest_path(1, 7))
    print(gr2.distance(1, 8))
    print(gr2.shortest_path(1, 8))
    print(gr2.distance(6, 1))
    print(gr2.shortest_path(6, 1))

    # print(gr2.reachable_with_dist(1))
    print()
    print(gr.mean_degree())
    print(gr.prob_degree())
    print(gr.mean_distances())
    print()
    print(gr.clustering_coef(1))
    print(gr.clustering_coef(2))

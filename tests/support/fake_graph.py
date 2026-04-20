class FakeGraph:
    def __init__(self, nodes, edges):
        self._nodes = list(nodes)
        self._edges = list(edges)

    def nodes(self):
        return list(self._nodes)

    def edges(self):
        return list(self._edges)

    def in_degree(self, node):
        return sum(1 for _, target in self._edges if target == node)

    def out_degree(self, node):
        return sum(1 for source, _ in self._edges if source == node)

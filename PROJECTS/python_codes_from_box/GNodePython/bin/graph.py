class Graph(object):

    def __init__(self, label):
        self.label = label
        self.nodes = []

    def addnode(self, label):
        print "got here"
        _counter = len(self.nodes)
        print "here 2" 
        print _counter
        self.nodes[_counter] = Node(label)


class Node(object):

    def __init__(self, label):
        self.label = label


class Edge(object):

    def __init__(self, node1, node2, directed = False, weight = 1.0):
        self.directed = directed
        self.weight = weight
        # here we could test the nodes exist
        self.nodes = (node1, node2)

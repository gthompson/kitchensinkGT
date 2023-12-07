import unittest
import graph
#import exceptions


class TestGraph(unittest.TestCase):
    
    def testNewGraph(self):
         
        #self.assertRaises(Exception, Graph)
        g = graph.Graph("graph1")
        n = graph.Node("node1")
        e = graph.Edge(1, 2)
        g.addnode("node2")
        self.assertTrue(len(g.nodes) == 1)   
    
        
if __name__ == '__main__':
    unittest.main()
    

        
        

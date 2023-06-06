import unittest
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import time

import eppstein
from epsnet.graph import generate_random_graph
from epsnet.utils import draw_graph

# key for sorting the list of (list, cost) tuples
def custom_key(item):
    l, n = item
    return (n, len(l), l)

def create_pos_weighted_graph():
    graph = nx.DiGraph() # Directed Graph
    graph.add_node('s', name = "source", index= 's')
    graph.add_node('t', name = "destination",index='t' )
    for i in range(3):
        for j in range(4):
                graph.add_node((i,j), index=(i,j), name = (i,j))
    edges = []
    edges.append(('s',(0,0),0))
    edges.append(((2,3),'t',0))
    
    edges.append(((0,0),(0,1),2))
    edges.append(((0,0),(1,0),13))
    
    edges.append(((0,1),(0,2),20))
    edges.append(((0,1),(1,1),27))   

    edges.append(((0,2),(0,3),14))        
    edges.append(((0,2),(1,2),14))            

    edges.append(((0,3),(1,3),15))                        
    
    edges.append(((1,0),(1,1),9))
    edges.append(((1,0),(2,0),15))                        
    
    edges.append(((1,1),(1,2),10))                        
    edges.append(((1,1),(2,1),20))
    
    edges.append(((1,2),(1,3),25))                        
    edges.append(((1,2),(2,2),12))
    
    edges.append(((1,3),(2,3),7))
    edges.append(((2,0), (2,1),18))
    edges.append(((2,1), (2,2),8))
    edges.append(((2,2), (2,3),11))
                                
                                
                                    
    graph.add_weighted_edges_from(edges)
    return graph

def create_neg_weighted_graph():
    graph = nx.DiGraph() # Directed Graph
    graph.add_node('s', name = "source", index= 's')
    graph.add_node('t', name = "destination",index='t' )
    for i in range(3):
        for j in range(4):
                graph.add_node((i,j), index=(i,j), name = (i,j))
    edges = []
    edges.append(('s',(0,0),0))
    edges.append(((2,3),'t',0))
    
    edges.append(((0,0),(0,1),-2))
    edges.append(((0,0),(1,0),-13))
    
    edges.append(((0,1),(0,2),-20))
    edges.append(((0,1),(1,1),-27))   

    edges.append(((0,2),(0,3),-14))        
    edges.append(((0,2),(1,2),-14))            

    edges.append(((0,3),(1,3),-15))                        
    
    edges.append(((1,0),(1,1),-9))
    edges.append(((1,0),(2,0),-15))                        
    
    edges.append(((1,1),(1,2),-10))                        
    edges.append(((1,1),(2,1),-20))
    
    edges.append(((1,2),(1,3),-25))                        
    edges.append(((1,2),(2,2),-12))
    
    edges.append(((1,3),(2,3),-7))
    edges.append(((2,0), (2,1),-18))
    edges.append(((2,1), (2,2),-8))
    edges.append(((2,2), (2,3),-11))
                                    
    graph.add_weighted_edges_from(edges)
    return graph

class TestKShortest(unittest.TestCase):

    def setUp(self):
        self.graph = nx.DiGraph()
        self.graph.add_nodes_from(["A","B","C","D","E","F","G","H","I"])
        self.graph.add_weighted_edges_from([
            ("A","B",2),
            ("A","C",0),
            ("B","D",2),
            ("B","E",6),
            ("C","E",8),
            ("C","F",3),
            ("D","G",5),
            ("E","G",3),
            ("E","H",4),
            ("F","H",4),
            ("G","I",5),
            ("H","I",6),
            ])
    
    # edited
    # python -m unittest
    def test_rnd_graph_eppstein(self):
        n = 10
        p = 0.5
        k = 10 # Number of paths to consider

        # lines = open("random/g_rand2", "r")
        # g_rand = nx.parse_edgelist(lines, create_using=nx.DiGraph(), nodetype=int)
        # lines.close()
        g_rand = generate_random_graph(n, p)
        # nx.write_edgelist(g_rand, 'g_rand2')
        # print(g_rand.nodes)
        # print(nx.is_directed_acyclic_graph(g_rand))
        # draw_graph(g_rand)
        # plt.show()


        # Get all available paths
        pathvec = list(nx.all_simple_paths(g_rand, 0, n-1))
        weightvec = [(x, nx.path_weight(g_rand, x, weight="weight")) for x in pathvec]
        pathvec_full = sorted(weightvec, key=lambda x:(x[1],tuple(x[0])))
        
        pathvec_eppstein = eppstein.k_shortest_paths(g_rand,0,n-1,k)

        pathvec_full.sort(key=custom_key)
        pathvec_eppstein.sort(key=custom_key)

        pathvec_full_cut = pathvec_full[:k]
        pathvec_eppstein_cut = pathvec_eppstein[:k]
        
        print()
        print(pathvec_full_cut)
        print(pathvec_eppstein_cut)

        self.assertEqual(pathvec_eppstein_cut,pathvec_full_cut)

    # edited
    def test_eppstein_waterman_example(self):
        val = eppstein.k_shortest_paths(self.graph,"A","I",2)
        ref = [(['A', 'C', 'F', 'H', 'I'], 13), (['A', 'B', 'D', 'G', 'I'], 14)]
        self.assertEqual(len(val),2)
        self.assertEqual(ref,val)

    # edited
    # def test_eppstein_neg_weight(self):
    #     graph = create_neg_weighted_graph()
    #     pathvec_eppstein =  eppstein.k_shortest_paths(graph,"s","t",3)

    #     pathvec = list(nx.all_simple_paths(graph, "s", "t"))
    #     weightvec = [(x, nx.path_weight(graph, x, weight="weight"))
    #                  for x in pathvec]
    #     pathvec_full = sorted(weightvec, key=lambda x:(x[1],tuple(x[0])))

    #     pathvec_full_cut = pathvec_full[:3]
    #     print(pathvec_full_cut)
    #     print(pathvec_eppstein)
    #     self.assertEqual(pathvec_eppstein,pathvec_full_cut)

    # edited
    def test_eppstein_pos_weight(self):
        graph = create_pos_weighted_graph()
        pathvec_eppstein =  eppstein.k_shortest_paths(graph,"s","t",10)

        pathvec = list(nx.all_simple_paths(graph, "s", "t"))
        weightvec = [(x, nx.path_weight(graph, x, weight="weight"))
                     for x in pathvec]
        pathvec_full = sorted(weightvec, key=lambda x:(x[1],tuple(x[0])))

        pathvec_full_cut = pathvec_full[:10]
        self.assertEqual(pathvec_eppstein,pathvec_full_cut)
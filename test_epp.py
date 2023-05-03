import unittest
import heapq
import networkx as nx
import matplotlib.pyplot as plt

# files in same dir
import eppstein
import test_graph
from epsnet.utils import draw_graph


class TestEpp(unittest.TestCase):

    # Set up and tear down
    @classmethod
    def setUpClass(cls):
        pass
    
    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        self.G = nx.DiGraph()
        self.G.add_weighted_edges_from([
            (2, 't', 7),
            (3, 't', 9),
            (3, 2, 10),
            ('s', 2, 15),
            ('s', 3, 11),
            ('s', 5, 6),
            (5, 6, 9),
            (6, 't', 14),
            (6, 3, 2),
        ])
        self.Gpred = {'t': [], 2: ['t'], 3: ['t'], 6: [3], 's': [3], 5: [6]}
        self.Gdist = {'t': 0, 2: 7, 3: 9, 6: 11, 's': 20, 5: 20}

        self.Gsidetrack = self.G.copy()
        self.Gsidetrack.add_edges_from([
            (2, 't', {'sidetrackCost': 0}),
            (3, 't', {'sidetrackCost': 0}),
            (3, 2, {'sidetrackCost': 8}),
            ('s', 2, {'sidetrackCost': 2}),
            ('s', 3, {'sidetrackCost': 0}),
            ('s', 5, {'sidetrackCost': 6}),
            (5, 6, {'sidetrackCost': 0}),
            (6, 't', {'sidetrackCost': 3}),
            (6, 3, {'sidetrackCost': 0}),
        ])

        self.GTree = nx.DiGraph()
        self.GTree.add_edges_from([
            (eppstein.STCedge(None, 's', {'weight': 0, 'sidetrackCost': 0}),
             eppstein.STCedge(3, 2, {'weight': 10, 'sidetrackCost': 8})),
            (eppstein.STCedge(None, 's', {'weight': 0, 'sidetrackCost': 0}),
             eppstein.STCedge('s', 2, {'weight': 15, 'sidetrackCost': 2})),
            (eppstein.STCedge(None, 's', {'weight': 0, 'sidetrackCost': 0}),
             eppstein.STCedge('s', 5, {'weight': 6, 'sidetrackCost': 6})),
            (eppstein.STCedge('s', 5, {'weight': 6, 'sidetrackCost': 6}),
             eppstein.STCedge(6, 't', {'weight': 14, 'sidetrackCost': 3})),
            (eppstein.STCedge('s', 5, {'weight': 6, 'sidetrackCost': 6}),
             eppstein.STCedge(3, 2, {'weight': 10, 'sidetrackCost': 8})),
        ])

        self.GHout = self.Gsidetrack.copy()
        # create the Hout elements
        vt = eppstein.Hout(self.Gsidetrack, 't')
        vt.root = None
        vt.heap = []
        v2 = eppstein.Hout(self.Gsidetrack, 2)
        v2.root = None
        v2.heap = []
        v3 = eppstein.Hout(self.Gsidetrack, 3)
        v3.root = eppstein.STCedge(3, 2, self.Gsidetrack.adj[3][2])
        v3.heap = []
        vs = eppstein.Hout(self.Gsidetrack, 's')
        vs.root = eppstein.STCedge('s', 2, self.Gsidetrack.adj['s'][2])
        vs.heap = [eppstein.STCedge('s', 5, self.Gsidetrack.adj['s'][5])]
        v5 = eppstein.Hout(self.Gsidetrack, 5)
        v5.root = None
        v5.heap = []
        v6 = eppstein.Hout(self.Gsidetrack, 6)
        v6.root = eppstein.STCedge(6, 't', self.Gsidetrack.adj[6]['t'])
        v6.heap = []
        self.GHout.nodes['t']['Hout'] = vt
        self.GHout.nodes[2]['Hout'] = v2
        self.GHout.nodes[3]['Hout'] = v3
        self.GHout.nodes['s']['Hout'] = vs
        self.GHout.nodes[5]['Hout'] = v5
        self.GHout.nodes[6]['Hout'] = v6
    
        # self.GH_G = self.GHout.copy()


        # Second test graph
        self.F = nx.DiGraph()
        self.F.add_weighted_edges_from([(0, 2, 5), (1, 2, 3), (0, 1, 2), (1, 3, 3), (3, 2, 1)])

        self.Fpred = {2: [], 0: [2, 1], 1: [2], 3: [2]}
        self.Fdist = {2: 0, 3: 1, 1: 3, 0: 5}

        self.Fsidetrack = nx.DiGraph()
        self.Fsidetrack.add_edges_from([
            (0, 2, {'weight': 5, 'sidetrackCost': 0}),
            (0, 1, {'weight': 2, 'sidetrackCost': 0}),
            (1, 2, {'weight': 3, 'sidetrackCost': 0}),
            (1, 3, {'weight': 3, 'sidetrackCost': 1}),
            (3, 2, {'weight': 1, 'sidetrackCost': 0}),
        ])

        self.FTree = nx.DiGraph()
        self.FTree.add_edges_from([
            (eppstein.STCedge(None, 0, {'weight': 0, 'sidetrackCost': 0}),
             eppstein.STCedge(1, 3, {'weight': 3, 'sidetrackCost': 1})),
        ])

        self.FHout = self.Fsidetrack.copy()
        v0 = eppstein.Hout(self.Fsidetrack, 0)
        v0.root = None
        v0.heap = []
        v1 = eppstein.Hout(self.Fsidetrack, 1)
        v1.root = eppstein.STCedge(1, 3, self.Fsidetrack.adj[1][3])
        v1.heap = []
        v2 = eppstein.Hout(self.Fsidetrack, 2)
        v2.root = None
        v2.heap = []
        v3 = eppstein.Hout(self.Fsidetrack, 3)
        v3.root = None
        v3.heap = []
        self.FHout.nodes[0]['Hout'] = v0
        self.FHout.nodes[1]['Hout'] = v1
        self.FHout.nodes[2]['Hout'] = v2
        self.FHout.nodes[3]['Hout'] = v3


        # From test_graph.py from the epsnet repo
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

    def tearDown(self):
        pass

    # Actual tests
    def test_calc_sidetrack_cost(self):
        val = self.G
        eppstein.calc_sidetrack_cost(val, self.Gdist)
        ref = self.Gsidetrack
        self.assertEqual(ref.adj,val.adj)

        val = self.F
        eppstein.calc_sidetrack_cost(val, self.Fdist)
        ref = self.Fsidetrack
        self.assertEqual(val.adj, ref.adj)


    def test_sidetrackEdge_path_tree(self):
        val = eppstein.sidetrackEdge_path_tree(self.Gsidetrack, self.Gpred, 's')
        ref = self.GTree
        self.assertEqual(ref.adj,val.adj)

        val = eppstein.sidetrackEdge_path_tree(self.Fsidetrack, self.Fpred, 0)

        ref = self.FTree
        self.assertEqual(val.adj, ref.adj)
    

    def test_Hout(self):
        for v in self.Gsidetrack.nodes():
            self.Gsidetrack.nodes[v]['Hout'] = eppstein.Hout(self.Gsidetrack,v)
        self.assertEqual(self.Gsidetrack.adj, self.GHout.adj)

        val = self.Fsidetrack
        for v in val.nodes():
            val.nodes[v]['Hout'] = eppstein.Hout(self.Fsidetrack,v)  
        self.assertEqual(val.adj, self.FHout.adj)


    def test_calc_H_G(self):
        # test explanation example
        val = eppstein.calc_H_G(self.GHout, self.Gpred, 't')

        v5 = val.nodes[5]['H_G']
        self.assertEqual(3,v5[0].root.strc)
        self.assertEqual(8,v5[1].root.strc)
        self.assertEqual([],v5[0].heap)
        self.assertEqual([],v5[1].heap)

        vS = val.nodes['s']['H_G']
        self.assertEqual(2,vS[0].root.strc)
        self.assertEqual(8,vS[1].root.strc)
        self.assertEqual(6,vS[0].heap[0].strc)
        self.assertEqual(1, len(vS[0].heap))
        self.assertEqual([],vS[1].heap)

        # test multiple shortest paths
        val = eppstein.calc_H_G(self.FHout, self.Fpred, 2)
        v0 = val.nodes[0]['H_G']
        v2 = val.nodes[2]['H_G']
        self.assertEqual(1,v0[0].root.strc)
        self.assertEqual([],v0[0].heap)
        self.assertEqual(v2,[])
        # for v in val.nodes():
        #     print(v, ":", val.nodes[v]['H_G'])


    # def test_P_to_Heap(self):
    # #     G = eppstein.calc_H_G(self.GHout, self.Gpred, 't')
    # #     H_G_dict = {}
    # #     for p in self.GTree.nodes: # p is of the form STCedge()
    # #         H_G_dict[p] = G.nodes[p.head]['H_G']
        
    # #     P = nx.DiGraph()
    # #     Proot = eppstein.prepare_and_augmentP(P, H_G_dict, 's')

    # #     print(next(iter(P.adj[Proot])))

    # #     H = nx.DiGraph()
    # #     Hroot = eppstein.P_to_Heap(H, P, Proot)


    # #     print()
    # #     for n in H.adj:
    # #         print(n, "->", H.adj[n])
    # #         print()

    # #     nx.draw(H)
    # #     plt.draw()
    # #     plt.show()
    #     G = eppstein.calc_H_G(self.FHout, self.Fpred, 2)
    #     H_G_dict = {}
    #     for p in self.FTree.nodes: # p is of the form STCedge()
    #         H_G_dict[p] = G.nodes[p.head]['H_G']

    #     P = nx.DiGraph()
    #     Proot = eppstein.prepare_and_augmentP(P, H_G_dict, 0)

    #     print(next(iter(P.adj[Proot])))

    #     H = nx.DiGraph()
    #     Hroot = eppstein.P_to_Heap(H, P, Proot)


    #     print()
    #     for n in H.adj:
    #         print(n, "->", H.adj[n])
    #         print()

    #     draw_graph(H)
    #     plt.show()


    def test_shortest_paths(self):
        # test explanation example
        sp = eppstein.k_shortest_paths(self.G, 's', 't', 6)
        ref = [
            (['s', 3, 't'], 20),
            (['s', 2, 't'], 22),
            (['s', 5, 6, 3, 't'], 26),
            (['s', 3, 2, 't'], 28),
            (['s', 5, 6, 't'], 29),
            (['s', 5, 6, 3, 2, 't'], 34)
        ]
        self.assertEqual(sp,ref)

        # test graph where not every node is reachable from source
        g = nx.DiGraph()
        g.add_weighted_edges_from([(0, 2, 5), (1, 2, 3)])
        ref = [([0, 2], 5)]
        val = eppstein.k_shortest_paths(g, 0, 2, 1)
        self.assertEqual(ref,val)

        # g.clear()
        # g.add_weighted_edges_from([(0, 2, 5), (1, 2, 3), (0, 1, 2)])
        g = self.F
        draw_graph(g)
        plt.show()
        print(eppstein.k_shortest_paths(g, 0, 2, 2))


if __name__ == '__main__':
    unittest.main()
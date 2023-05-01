import unittest
import heapq
import networkx as nx
import matplotlib.pyplot as plt

# files in same dir
import eppstein
import test_graph


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
    
        self.GH_G = self.GHout.copy()


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
        val = eppstein.calc_sidetrack_cost(self.G, self.Gdist)
        ref = self.Gsidetrack
        self.assertEqual(ref.adj,val.adj)

    def test_sidetrackEdge_path_tree(self):
        val = eppstein.sidetrackEdge_path_tree(self.Gsidetrack, self.Gpred, 's')
        ref = self.GTree
        self.assertEqual(ref.adj,val.adj)
    
    def test_Hout(self):
        for v in self.Gsidetrack.nodes():
            self.Gsidetrack.nodes[v]['Hout'] = eppstein.Hout(self.Gsidetrack,v)
        
        # for v in self.Gsidetrack.nodes():
        #     print(v, ":", self.Gsidetrack.nodes[v]['Hout'])
        # print()
        # for v in self.GHout.nodes():
        #     print(v, ":", self.GHout.nodes[v]['Hout'])

        self.assertEqual(self.Gsidetrack.adj, self.GHout.adj)

    def test_calc_H_G(self):
        val = eppstein.calc_H_G(self.GHout, self.Gpred, 't')

        # for v in val.nodes():
        #     print(v, ":", val.nodes[v]['H_G'])


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

    # def test_P_to_Heap(self):
    #     G = eppstein.calc_H_G(self.GHout, self.Gpred, 't')
    #     H_G_dict = {}
    #     for p in self.GTree.nodes: # p is of the form STCedge()
    #         H_G_dict[p] = G.nodes[p.head]['H_G']
        
    #     P = nx.DiGraph()
    #     Proot = eppstein.prepare_and_augmentP(P, H_G_dict, 's')

    #     print(next(iter(P.adj[Proot])))

    #     H = nx.DiGraph()
    #     Hroot = eppstein.P_to_Heap(H, P, Proot)


    #     print()
    #     for n in H.adj:
    #         print(n, "->", H.adj[n])
    #         print()

    #     nx.draw(H)
    #     plt.draw()
    #     plt.show()

    def test_shortest_paths(self):
        sp = eppstein.shortest_paths(self.G, 's', 't', 6)
        ref = [
            (['s', 3, 't'], 20),
            (['s', 2, 't'], 22),
            (['s', 5, 6, 3, 't'], 26),
            (['s', 3, 2, 't'], 28),
            (['s', 5, 6, 't'], 29),
            (['s', 5, 6, 3, 2, 't'], 34)
        ]
        self.assertEqual(sp,ref)

if __name__ == '__main__':
    unittest.main()
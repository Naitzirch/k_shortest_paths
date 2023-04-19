import unittest
import networkx as nx
import eppstein

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

        self.Gsidetrack = self.G
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
            ((), (3, 2, 8)),
            ((), ('s', 2, 2)),
            ((), ('s', 5, 6)),
            (('s', 5, 6), (6, 't', 3)),
            (('s', 5, 6), (3, 2, 8)),
        ])
    
    def tearDown(self):
        pass

    # Actual tests
    def test_calc_sidetrack_cost(self):
        val = eppstein.calc_sidetrack_cost(self.G, self.Gdist)
        ref = self.Gsidetrack
        self.assertEqual(ref,val)

    def test_sidetrackEdge_path_tree(self):
        val = eppstein.sidetrackEdge_path_tree(self.Gsidetrack, self.Gpred, 's')
        ref = self.GTree
        self.assertEqual(ref.adj,val.adj)

if __name__ == '__main__':
    unittest.main()
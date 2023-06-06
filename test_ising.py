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

import ast
def parse_tuple(string):
    try:
        s = ast.literal_eval(str(string))
        if type(s) == tuple:
            return s
        return
    except:
        return

def set_subset(dest, n):
    # Distribute subsets according to the sites (this makes drawing easier)
    nodes = []
    for node in dest.nodes():
        if node != 's' and node != 't':
            nodes.append(parse_tuple(node))
    
    subsetdict = {str((x,i)):{"subset":x} for (x,i) in nodes}
    nx.set_node_attributes(dest,subsetdict)
    dest.nodes['s']['subset'] = -1
    dest.nodes['t']['subset'] = n

def draw_system_graph(graph,ax=None):
    """Draw the connectivity graph for the nsite system.
    This graph is fully connected between the different sites.
    The number of nodes in each layer is given by the number of points in the respective eps-net on the site.

    Args:
        ax (plt.axis, optional): Matplotlib axis to use for drawing. Defaults to None.
    """
    f, axvec = plt.subplots(1,1)
    pos = nx.multipartite_layout(graph)
    labels = nx.get_edge_attributes(graph,'weight')
    new_labels = {}
    for k, v in labels.items():
        new_labels[k] = round(v, 2)

    pred, dist = nx.dijkstra_predecessor_and_distance(graph.reverse(copy=True), 't')
    #get edges in pred
    cedges = []
    for u in pred:
        for v in pred[u]:
            cedges.append((u, v))

    edge_color = []
    for e in graph.edges():
        if e in cedges:
            edge_color.append('red')
        else:
            edge_color.append('black')

    nx.draw_networkx(graph, pos, ax=ax,with_labels=True)
    nx.draw_networkx_edges(graph, pos, width=2.0, edge_color=edge_color)
    # nx.draw_networkx_edge_labels(graph, pos, ax=ax,edge_labels=new_labels)

class TestKShortest(unittest.TestCase):

    def test_eppstein_ising_graphs(self):
        k = 10 # Number of paths to consider

        #g = nx.read_graphml("epsnet/IsingModel/ising_nsites_2_npoints_4.gml")

        #g = nx.read_graphml("epsnet/IsingModel/ising_nsites_3_npoints_4.gml")
        #g = nx.read_graphml("epsnet/IsingModel/ising_nsites_4_npoints_4.gml")
        #g = nx.read_graphml("epsnet/IsingModel/ising_nsites_5_npoints_4.gml")
        #g = nx.read_graphml("epsnet/IsingModel/ising_nsites_5_npoints_6.gml")
        #g = nx.read_graphml("epsnet/IsingModel/ising_nsites_5_npoints_10.gml")
        #g = nx.read_graphml("epsnet/IsingModel/ising_nsites_5_npoints_20.gml")

        
        # perturbated
        #g = nx.read_graphml("epsnet/IsingModel/perturbated/p_ising_nsites_2_npoints_4.gml")
        #g = nx.read_graphml("epsnet/IsingModel/perturbated/p_ising_nsites_3_npoints_4.gml")
        #g = nx.read_graphml("epsnet/IsingModel/perturbated/p_ising_nsites_4_npoints_4.gml")
        #g = nx.read_graphml("epsnet/IsingModel/perturbated/p_ising_nsites_5_npoints_4.gml")
        g = nx.read_graphml("epsnet/IsingModel/perturbated/p_ising_nsites_5_npoints_6.gml")
        #g = nx.read_graphml("epsnet/IsingModel/perturbated/p_ising_nsites_5_npoints_10.gml")
        #g = nx.read_graphml("epsnet/IsingModel/perturbated/p_ising_nsites_5_npoints_20.gml")
        
        #set_subset(g, 10)
        #draw_system_graph(g)
        #plt.show()

        ##############################
        pathvec = list(nx.all_simple_paths(g, 's', 't'))
        weightvec = [(x, nx.path_weight(g, x, weight="weight")) for x in pathvec]
        pathvec_full = sorted(weightvec, key=lambda x:(x[1],tuple(x[0])))
    
        pathvec_full.sort(key=custom_key)
        pathvec_full = [(x, round(y, 7)) for x, y in pathvec_full]
        pathvec_full_cut = pathvec_full[:k]

        with open(r'exclude/output/brutefor_n5.txt', 'w') as fp:
            for line in pathvec_full_cut:
                fp.write(f"{line}\n")
        ##############################

        #############################
        pathvec_eppstein = eppstein.k_shortest_paths(g,'s','t',k)

        pathvec_eppstein.sort(key=custom_key)
        pathvec_eppstein = [(x, round(y, 7)) for x, y in pathvec_eppstein]
        pathvec_eppstein_cut = pathvec_eppstein[:k]

        with open(r'exclude/output/eppstein_n5.txt', 'w') as fp:
            for line in pathvec_eppstein_cut:
                fp.write(f"{line}\n")
        ##############################

        print()
        print(pathvec_full)
        print(pathvec_eppstein_cut)

        self.assertEqual(pathvec_eppstein_cut,pathvec_full_cut)
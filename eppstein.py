''' Notes
path lenghts must be numerical
'''

''' Self Notes

dijkstra_predecessor_and_distance:
    pred and dist are dictionaries:
    sorted(pred.items())
    >>> [(0, []), (1, [0]), (2, [1]), (3, [2]), (4, [3])]
    sorted(dist.items())
    >>> [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]
    The list of predecessors contains more than one element only when there are more than one shortest paths to the key node.


'''


import math
import copy
import heapq
import networkx as nx



# # get_sidetrack_edges routines will run infinite if G has cycles... -> solution with dist[] dictionary and counter? or seen unseen counter?
# # DFS
# def get_sidetrack_edges_DFS(G, pred, src):
#     v = src
#     temp_dict = {}
#     # iterate over all nodes on the shortest path from src to dst
#     while v != None:
#         # iterate over all out edges of v
#         for nbr in G.adj[v]:
#             stc = G.adj[v][nbr]['sidetrackCost']
#             if stc != 0:
#                 temp_dict[(v, nbr, stc)] = get_sidetrack_edges_DFS(G, pred, nbr)
#         v = pred.get(v)

#     return temp_dict

# # BFS
# def get_sidetrack_edges_BFS(G, pred, src):
#     v = src
#     temp_dict = {}
#     # iterate over all nodes on the shortest path from src to dst
#     while v != None:
#         # iterate over all out edges of v
#         for nbr in G.adj[v]:
#             stc = G.adj[v][nbr]['sidetrackCost']
#             if stc != 0:
#                 temp_dict[(v, nbr, stc)] = None
#         v = pred.get(v)

#     for x in temp_dict:
#         temp_dict[x] = get_sidetrack_edges_BFS(G, pred, nbr)

#     return temp_dict

class STCedge:
    def __init__(self, tail, head, attr):
        self.tail = tail
        self.head = head
        self.attr = attr
        self.strc = attr['sidetrackCost']
    def __lt__(self, other):
        return self.strc < other.strc
    def __str__(self):
        return f"{self.strc}"

# Heap of the out edges of v
class Hout:
    def __init__(self, G, v):
        self.smallest = math.inf
        self.root = None            # outgoing edge with smallest stc value
        self.heap = []
        for u in G.adj[v]:
            e = STCedge(v, u, G.adj[v][u])

            # Hout consist only of non-0 stc edges
            if e.strc == 0:
                continue

            if e.strc < self.smallest:
                # add current root to internal heap unless root is None
                if self.root != None:
                    heapq.heappush(self.heap, self.root)
                
                # update root
                self.smallest = e.strc
                self.root = e
            else:
                heapq.heappush(self.heap, e)
    def __lt__(self, other):
        return self.root < other.root
    def __str__(self):
        return f"{self.root}, {self.heap}"


# get the sidetrack edges with tails on the shortest path from u to dst
# G: the original graph
# pred: next node on the shortest path tree from u to dst
# u: source node
# STree: The Sidetrack Sequence path tree to be constructed
# prevNode: the previous node in STree to which the new found sidetrack edges should be appended
def get_sidetrack_edges_DFS(G, pred, u, STree, prevNode):
    v = [u]
    # iterate over all nodes on the shortest path from src to dst
    while v != []:
        v = v[0]
        for nbr in G.adj[v]:
            stc = G.adj[v][nbr]['sidetrackCost']
            if stc != 0:
                STree.add_edge(prevNode, (v, nbr, stc))
                STree = get_sidetrack_edges_DFS(G, pred, nbr, STree, (v, nbr, stc))
        v = pred.get(v)
    
    return STree

# get the sidetrack edges with tails on the shortest path from u to dst
# G: the original graph
# pred: next node on the shortest path tree from u to dst
# u: source node
# STree: The Sidetrack Sequence path tree to be constructed
# prevNode: the previous node in STree to which the new found sidetrack edges should be appended
def get_sidetrack_edges_BFS(G, pred, u, STree, prevNode):
    v = [u]
    # iterate over all nodes on the shortest path from src to dst
    while v != []:
        v = v[0]
        for nbr in G.adj[v]:
            stc = G.adj[v][nbr]['sidetrackCost']
            if stc != 0:
                STree.add_edge(prevNode, (v, nbr, stc))
        v = pred.get(v)
    
    for child in STree.adj[prevNode]:
        STree = get_sidetrack_edges_BFS(G, pred, child[1], STree, child)

    return STree

# nodes in the sidetrack edge path tree will be of the form
# (tail, head, sidetrackCost)
# the path from root to node defines the sequence of sidetrack edges
# i.e. each node is lastsidetrack(p)
def sidetrackEdge_path_tree(G, pred, src, dst):
    # create a tree of sequences of sidetrack edges, denoting paths in G
    STree = nx.DiGraph()

    # add the empty sequence as root
    STree.add_node((None, dst, 0))

    return get_sidetrack_edges_DFS(G, pred, src, STree, ())


def calc_sidetrack_cost(G, dist):
    # iterate over all nodes in adjecency list form ( should take O(m) )
    for u, nbrs in G.adj.items():
        # remove nodes from which dst can't be reached
        if dist.get(u) == None:
            G.remove_node(u)
        
        # Add sidetrack costs as attribute of edges to the edges in G
        # [and connections from u to head(sidetrack-edge) for each side-track edge on the shortest path from u to t]?
        # delta(e) = l(e) + dist(head(e), t) - dist(tail(e),t)
        else:
            for nbr, eattr in nbrs.items():
                G[u][nbr]['sidetrackCost'] = eattr['weight'] + dist[nbr] - dist[u]
    return G

# DFS to create H_G heaps for each vertex v
# ordered by the value of the roots of each Hout heap on the shortest path from v to t
def calc_H_G_next(R, pred, prevNode):
    for v in R.adj[prevNode]:
        if prevNode in pred[v]:
            Hout_v = R.nodes[v]['Hout']
            h = R.nodes[v]['H_G'] = copy.copy(R.nodes[prevNode]['H_G'])
            if Hout_v.root != None:
                heapq.heappush(h, Hout_v)
            calc_H_G_next(R, pred, v)
                
# For each vertex v in G, creates a heap H_G of all Hout heaps on the path from v to t
# ordered by value of the roots of each Hout heap
def calc_H_G(G, pred, dst):
    R = G.reverse(copy=True)
    h = R.nodes[dst]['H_G'] = []
    Hout_dst = R.nodes[dst]['Hout']
    
    if Hout_dst.root != None:
        heapq.heappush(h, Hout_dst)
    
    calc_H_G_next(R, pred, dst)

    return R.reverse(copy=True)

# DON'T FORGET INDEX OUT OF RANGE ERRORS
def Hout_DFS(P, h, i):
    P.add_edge(h[i], h[2*i+1])
    Hout_DFS(P, h, i)
    P.add_edge(h[i], h[2*i+2])
    Hout_DFS(P, h, i)

def HoutHeap_DFS(P, h, i):
    # Add the two edges leading to other Hout heaps
    P.add_edge(h[i].root, h[2*i+1].root)
    HoutHeap_DFS(P, h, 2*i+1)
    P.add_edge(h[i].root, h[2*i+2].root)
    HoutHeap_DFS(P, h, 2*i+2)

    # Add its own heap
    P.add_edge(h[i].root, h[i].heap[0])
    Hout_DFS(P, h[i], 0)

# Transform all the heaps into nodes in 1 digraph
def prepareP(P, H_G_list):
    # iterate over H_G
    for HoutHeap in H_G_list:
        # HoutHeap = heap of Hout objects
        # Traverse the Hout heap, add edges from the root elements to the lower roots
        if HoutHeap != []:
            HoutHeap_DFS(P, HoutHeap, 0)

def augmentP(P):
    pass

# G:    a networkx DiGraph
# src:  the source node
# dst:  the destination node
# k:    the number of paths
def shortest_paths(G, src, dst, k):
    # Calculate single-destination shortest path tree; this is the
    # same as a single-source shortest path tree in G-reverse
    R = G.reverse(copy=True)

    # Essentially calculate min-weight spanning tree but just remember
    # predecessor of node on shortest path to dst and min distance to each node from
    # a given source
    pred, dist = nx.dijkstra_predecessor_and_distance(R, dst) # ( should take O(m + nlogn) )

    # Check if src and dst are connected
    if dist.get(src) == None:
        return []

    # Calculate sidetrack costs for every edge and add them as attribute to the edge
    G = calc_sidetrack_cost(G, dist)
    
    # Create a path tree with sidetrack(p) sequences S
    # where the parent of any path p is prefpath(p)
    # This tree will be heap-ordered by Lemma 3: l(p) >= l(prefpath(p))
    # (the empty sequence will be the root)
    pathTree = sidetrackEdge_path_tree(G, pred, src, dst)

    # Calculate Hout(v) for every vertex and add the resulting heap as attribute to v
    for v in G.nodes():
        G.nodes[v]['Hout'] = Hout(G,v)

    # Calculate H_G(v) for every vertex, a balanced heap of the roots of Hout on the path from v to t
    # where every root also points to the rest its original Hout heap
    G = calc_H_G(G, pred, dst)

    # Retrieve H_G(head(lastsidetrack(p))) for each node p in pathTree
    # such that we can build P(G)
    H_G_list = []
    for p in pathTree.nodes: # p is of the form (v, u, stc)
        H_G_list.append(G.nodes[p[1]]['H_G'])

    P = nx.DiGraph()
    prepareP(P, H_G_list) # results in a graph of STCedge elements

    augmentP(P)
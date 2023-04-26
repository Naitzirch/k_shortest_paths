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
    def __hash__(self):
        return hash((self.tail, self.head, self.attr['weight'], self.strc))
    def __eq__(self, other) : 
        return self.__dict__ == other.__dict__
    def __ne__(self, other):
        return not(self.__dict__ == other.__dict__)

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
            e = STCedge(v, nbr, G.adj[v][nbr])
            if e.strc != 0:
                STree.add_edge(prevNode, e)
                STree = get_sidetrack_edges_DFS(G, pred, nbr, STree, e)
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
def sidetrackEdge_path_tree(G, pred, src):
    # create a tree of sequences of sidetrack edges, denoting paths in G
    STree = nx.DiGraph()

    # add the empty sequence as root
    root = STCedge(None, src, {'weight': 0, 'sidetrackCost': 0})
    STree.add_node(root)

    return get_sidetrack_edges_DFS(G, pred, src, STree, root)


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
def Hout_DFS(P, h, i, H_G_dict):
    P.add_edge(h[i], h[2*i+1], weight=(h[2*i+1].strc - h[i].strc) ) # For edge (u, v) in D(G), add as edge weight: d(v) - d(u)
    Hout_DFS(P, h, i, H_G_dict)
    P.add_edge(h[i], h[2*i+2], weight=(h[2*i+2].strc - h[i].strc) )
    Hout_DFS(P, h, i, H_G_dict)

    # Add an edge from p=h[i] (that corresponds to (u, w)) to h(w) (aka h(p.head)) with weight d(h(p.head))
    # (= edge from p to H_G_dict[p][0].root)
    h_w = H_G_dict[h[i]][0].root
    P.add_edge(h[i], h_w, weight=h_w.strc)

def HoutHeap_DFS(P, h, i, H_G_dict):
    # Add the two edges leading to other Hout heaps
    p = h[i].root

    # Left Hout child
    c1 = h[2*i+1].root
    P.add_edge(p, c1, weight=(c1.strc - p.strc) ) # For edge (u, v) in D(G), add as edge weight: d(v) - d(u)
    HoutHeap_DFS(P, h, 2*i+1, H_G_dict)

    # Right Hout child
    c2 = h[2*i+2].root
    P.add_edge(p, c2, weight=(c2.strc - p.strc) )
    HoutHeap_DFS(P, h, 2*i+2, H_G_dict)

    # Add its own (STCedge) heap (inner child)
    P.add_edge(h[i].root, h[i].heap[0], weight=(h[i].heap[0].strc - h[i].root.strc))
    Hout_DFS(P, h[i].heap, 0, H_G_dict)

    # Add an edge from p (that corresponds to (u, w)) to h(w) (aka h(p.head)) with weight d(h(p.head))
    # (= edge from p to H_G_dict[p][0].root)
    h_w = H_G_dict[p][0].root
    P.add_edge(p, h_w, weight=h_w.strc)

# Transform all the heaps into nodes in 1 digraph
def prepare_and_augmentP(P, H_G_dict, src):

    # augmentation: Add a root node
    empty_seq = STCedge(None, src, {'weight': 0, 'sidetrackCost': 0})
    # and connect it to h(src)
    h_s = H_G_dict[empty_seq][0].root
    P.add_edge(empty_seq, h_s, weight=h_s.strc)

    # iterate over H_G
    for l in H_G_dict: # last sidetrack edge in sequence for path p
        # H_G = heap of Hout objects
        H_G = H_G_dict[l]
        # Traverse the Hout heap, add edges from the root elements to the lower roots
        if H_G != []:
            HoutHeap_DFS(P, H_G, 0, H_G_dict)
    
    return empty_seq


class EHeapElement:
    def __init__(self, seq, weight):
        self.seq = seq
        self.weight = weight
    def __lt__(self, other):
        return self.weight < other.weight
    def __hash__(self):
        return hash((self.seq, self.weight))
    def __eq__(self, other) : 
        return self.__dict__ == other.__dict__
    def __ne__(self, other):
        return not(self.__dict__ == other.__dict__)

# node: roote node to start BFS from
def P_to_Heap(H, P, node):

    n1 = next(iter(P.adj[node]))
    Er = EHeapElement([node], 0)
    En = EHeapElement([node, n1], P.adj[node][n1]['weight'])
    H.add_edge(Er, En)

    Er = EHeapElement([node], P.adj[node][n1]['weight'])

    # BFS P
    queue = []
    queue.append( (node, n1, Er) )

    while queue:
        l = queue.pop(0)
        l_t = l[0] # tail
        l_h = l[1] # head
        l_s = l[2].seq
        l_w = l[2].weight

        for child in P.adj[l_h]:
                if P.adj[l_h][child]['cross_edge'] == True:
                    m = EHeapElement(l_s.append(l_t), l_w)
                    n = EHeapElement(l_s.append(l_t).append(child), (l_w + P.adj[l_h][child]['weight']))
                    o = EHeapElement(l_s.append(l_t), (l_w + P.adj[l_h][child]['weight']))
                    H.add_edge(m, n)
                    queue.append( (l_h, child, o) )
                else:
                    m = EHeapElement(l_s.append(child), (l_w + P.adj[l_h][child]['weight']))
                    H.add_edge(l[2], m)
                    n = EHeapElement(l_s, (l_w + P.adj[l_h][child]['weight']))
                    queue.append( (l_h, child, n) )
    return Er

# H, a root, heap tuple
def pop_from_H(root, H):
    v = root
    last_root = v
    while H.adj[v]:
        # Get minimum of child of root
        m = min(H.adj[v])
        
        # remove edge between root and minimum of its childeren
        H.remove_edge(v, m)

        # save original childeren of root and minimum child of root
        r_childeren = copy.copy(H.adj[v])
        m_childeren = copy.copy(H.adj[m])

        # remove the root and the minimum child
        H.remove_node(v)
        H.remove_node(m)

        # add the original childeren of root as childeren of the minimum
        H.add_edges_from(list((m, w) for w in r_childeren))
        # add the original childeren of minimum as childeren of the root
        H.add_edges_from(list((v, n) for n in m_childeren))

        # add root as child of minimum
        H.add_edge(m, v)

        # add from last root to the minimum
        if last_root != root:
            H.add_edge(last_root, m)
        last_root = m

    H.remove_node(v)
    return root


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
    G = calc_H_G(G, pred, dst) # For every vertex H_G is added as an attribute

    # Retrieve H_G(head(lastsidetrack(p))) for each node p in pathTree
    # such that we can build P(G)
    H_G_dict = {}
    for p in pathTree.nodes: # p is of the form STCedge()
        H_G_dict[p] = G.nodes[p.head]['H_G']

    P = nx.DiGraph()
    Proot = prepare_and_augmentP(P, H_G_dict, src) # results in a graph of STCedge elements
    # P is now the completed path graph P(G), its root is returned by the function above

    # Lastly, P(G) will need to be transformed into a 4-heap H(G), so that the nodes
    # in H(G) represent paths in G. H(G) is constructed by forming a node for each path in
    # P(G) rooted at Proot.
    H = nx.DiGraph()
    Hroot = P_to_Heap(H, P, Proot)
    H = (Hroot, H)

    # To find the k shortest paths, pop k times from H and append to list
    paths = []
    for _ in range(k):
        paths.append(pop_from_H(*H))

    return paths
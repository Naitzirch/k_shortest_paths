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
from epsnet.utils import draw_graph
import matplotlib.pyplot as plt
import time


class STCedge:
    def __init__(self, tail, head, attr):
        self.tail = tail
        self.head = head
        self.attr = attr
        self.strc = attr['sidetrackCost']
        self.mark = 0
    def __lt__(self, other):
        return self.strc < other.strc
    def __str__(self):
        return f'STCedge(\'{self.tail}\', \'{self.head}\' {self.strc}, \'{self.mark}\')'
        #return f"{self.strc}"
    def __repr__(self):
        #return f"STCedge({self.strc})"
        return f'STCedge(\'{self.tail}\', \'{self.head}\' {self.strc}, \'{self.mark}\')'
    def __hash__(self):
        return hash((self.tail, self.head, self.attr['weight'], self.strc, self.mark))
    def __eq__(self, other):
        if isinstance(other, STCedge) and self.__class__ == other.__class__:
            return self.__dict__ == other.__dict__
        else:
            return False
    def __ne__(self, other):
        if not(isinstance(other, STCedge)): return True
        if not(self.__class__ == other.__class__): return True
        return not(self.__dict__ == other.__dict__)

# Heap of the out edges of v
class Hout:
    def __init__(self, G, v):
        self.smallest = math.inf
        self.root = None            # outgoing edge with smallest stc value
        self.heap = []
        for u in G.adj[v]:
            e = STCedge(v, u, G.adj[v][u])

            # Hout consists only of non-0 stc edges
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
    def __repr__(self):
        return f'Hout({self.root}, {self.heap})'


# # get_sidetrack_edges routines will run infinite if G has cycles... -> solution with dist[] dictionary and counter? or seen unseen counter?

# get the sidetrack edges with tails on the shortest path from u to dst
# G: the original graph
# pred: next node on the shortest path tree from u to dst
# u: source node
# STree: The Sidetrack Sequence path tree to be constructed
# prevNode: the previous node in STree to which the new found sidetrack edges should be appended
def get_sidetrack_edges_DFS(G, pred, u, STree, prevNode):
    vl = [u]
    # iterate over all nodes on all shortest paths from src to dst
    while vl != []:
        for v in vl:
            for nbr in G.adj[v]:
                e = STCedge(v, nbr, G.adj[v][nbr])
                if e.strc != 0:
                    STree.add_edge(prevNode, e)
                    STree = get_sidetrack_edges_DFS(G, pred, nbr, STree, e)
        vl = pred.get(v)
    
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

def clean_graph(G, dist):
    tbr = []
    for u in G.nodes:
        if dist.get(u) == None:
            tbr.append(u)
    G.remove_nodes_from(tbr)

def calc_sidetrack_cost(G, dist):
    # iterate over all nodes in adjecency list form ( should take O(m) )
    for u, nbrs in G.adj.items():
        # Add sidetrack costs as attribute of edges to the edges in G
        # [and connections from u to head(sidetrack-edge) for each side-track edge on the shortest path from u to t]?
        # delta(e) = l(e) + dist(head(e), t) - dist(tail(e),t)
        for nbr, eattr in nbrs.items():
            G[u][nbr]['sidetrackCost'] = eattr['weight'] + dist[nbr] - dist[u]


#######################
# keep track of how many times an edge has been copied and made unique
GLOB_edge_unique_copy_counter = 0

def make_unique(Hout):
    global GLOB_edge_unique_copy_counter
    GLOB_edge_unique_copy_counter += 1
    Hout.root.mark = GLOB_edge_unique_copy_counter
    return

'''H_G heap operations'''
'''Partially adapted from: https://gist.github.com/mumbleskates/0ef75bf3f25d0faeecc73ddb9373ea75'''

def parent_ix(ix):
    return (ix - 1) // 2

def children_ix(ix):
    return (ix * 2 + 1, ix * 2 + 2)

def mark_path_above(h, ix):
    while ix:   # update the rest of the path above inserted element
        p_ix = parent_ix(ix)
        make_unique(h[p_ix])
        ix = p_ix

def mark_path_below(h, p): # marks path below p including p itself
    make_unique(h[p])
    for c in children_ix(p):
        if c < len(h):
            mark_path_below(h, c)

def mark_path_below_init(h, p_ix, ix):
    c1, c2 = children_ix(p_ix)
    other_child = 0
    if ix == c1:
        other_child = c2
    else:
        other_child = c1
    
    if other_child < len(h):
        mark_path_below(h, other_child) # marks path below other_child including other_child itself

def sift_towards_root(h, ix):
    current = h[ix]
    while ix:   # stop when ix is 0, the root
        p_ix = parent_ix(ix)
        parent = h[p_ix]
        if min(parent, current) is parent:
            break
        make_unique(parent)
        # mark path below the other child of the parent because we're gonna swap them
        mark_path_below_init(h, p_ix, ix)
        h[p_ix], h[ix], ix = current, parent, p_ix
    mark_path_above(h, ix)

# push Hout element onto heap h
# mark all root elements that are updated by swap operations with a unique identifier
def heappush(h, Hout):
    if Hout.root == None:
        raise Exception("Hout root cannot be None")
    h.append(Hout)
    sift_towards_root(h, len(h)-1)

#######################

# BFS to create H_G heaps for each vertex v
# ordered by the value of the roots of each Hout heap on the shortest path from v to t
def calc_H_G_BFS(R, pred, dst):
    queue = [dst]
    while queue:
        m = queue.pop(0)
        for nbr in R.adj[m]:
            if m in pred[nbr]:
                h = copy.deepcopy(R.nodes[m]['H_G'])
                if R.nodes[nbr]['H_G'] == []:
                    # if H_G still empty: add H_G of previous node and Hout of nbr
                    Hout_nbr = R.nodes[nbr]['Hout']
                    if Hout_nbr.root != None:
                        heappush(h, Hout_nbr)
                    R.nodes[nbr]['H_G'] = h
                    queue.append(nbr)
                else:
                    # if H_G was already updated through another path, only add the additional
                    # H_G of the current path by merging the heaps
                    for Hout in h:
                        # will take O(N/2(log(N/2) + log(N))), can we do better? O(N/2)? using merge/heapify?
                        make_unique(Hout)
                        heappush(R.nodes[nbr]['H_G'], Hout)


# For each vertex v in G, creates a heap H_G of all Hout heaps on the path from v to t
# ordered by value of the roots of each Hout heap
def calc_H_G(G, pred, dst):
    R = G.reverse(copy=True)
    nx.set_node_attributes(R, [], 'H_G') # set H_G of each node to [] so we can later merge them
    h = R.nodes[dst]['H_G']
    Hout_dst = R.nodes[dst]['Hout']
    
    if Hout_dst.root != None:
        heappush(h, Hout_dst)
    
    calc_H_G_BFS(R, pred, dst)

    return R.reverse(copy=True)


def Hout_DFS(P, h, i, G):
    if 2*i+1 < len(h):
        P.add_edges_from([ (h[i], h[2*i+1], {'weight': (h[2*i+1].strc - h[i].strc), 'cross_edge': False}) ]) # For edge (u, v) in D(G), add as edge weight: d(v) - d(u)
        Hout_DFS(P, h, 2*i+1, G)
    if 2*i+2 < len(h):
        P.add_edges_from([ (h[i], h[2*i+2], {'weight': (h[2*i+2].strc - h[i].strc), 'cross_edge': False}) ])
        Hout_DFS(P, h, 2*i+2, G)

    # CROSS_EDGE
    # Add an edge from p=h[i] (that corresponds to (u, w)) to h(w) (aka h(p.head)) with weight d(h(p.head))
    # (= edge from p to H_G_dict[p][0].root)
    h_w = G.nodes[h[i].head]['H_G']
    if h_w != []:
        h_w = h_w[0].root
        P.add_edges_from([ (h[i], h_w, {'weight': h_w.strc, 'cross_edge': True}) ])

def HoutHeap_DFS(P, h, i, G):
    # Add the two edges leading to other Hout heaps
    p = h[i].root

    # Left Hout child
    if 2*i+1 < len(h):
        c1 = h[2*i+1].root
        P.add_edges_from([ (p, c1, {'weight': (c1.strc - p.strc), 'cross_edge': False}) ]) # For edge (u, v) in D(G), add as edge weight: d(v) - d(u)
        HoutHeap_DFS(P, h, 2*i+1, G)

    # Right Hout child
    if 2*i+2 < len(h):
        c2 = h[2*i+2].root
        P.add_edges_from([ (p, c2, {'weight': (c2.strc - p.strc), 'cross_edge': False}) ])
        HoutHeap_DFS(P, h, 2*i+2, G)

    # Add its own (STCedge) heap (inner child)
    if h[i].heap != []:
        P.add_edges_from([ (h[i].root, h[i].heap[0], {'weight': (h[i].heap[0].strc - h[i].root.strc), 'cross_edge': False}) ])
        Hout_DFS(P, h[i].heap, 0, G)

    # CROSS_EDGE
    # Add an edge from p (that corresponds to (u, w)) to h(w) (aka h(p.head)) with weight d(h(p.head))
    # (= edge from p to H_G_dict[p][0].root)
    h_w = G.nodes[p.head]['H_G']
    if h_w != []:
        h_w = h_w[0].root
        P.add_edges_from([ (p, h_w, {'weight': h_w.strc, 'cross_edge': True}) ])

# Transform all the heaps into nodes in 1 digraph
def prepare_and_augmentP(P, G, src):

    # augmentation: Add a root node
    empty_seq = STCedge(None, src, {'weight': 0, 'sidetrackCost': 0})
    # and connect it to h(src)
    h_s = G.nodes[src]['H_G']
    if h_s == []:
        P.add_node(empty_seq)
        return empty_seq
    else:
        h_sr = h_s[0].root
        P.add_edges_from([ (empty_seq, h_sr, {'weight': h_sr.strc, 'cross_edge': False}) ])

    # iterate over all nodes in G because all nodes are in the shortest path tree
    for l in G.nodes: # last sidetrack edge in sequence for path p
        # H_G = heap of Hout objects
        H_G = G.nodes[l]['H_G']
        # Traverse the Hout heap, add edges from the root elements to the lower roots
        if H_G != []:
            HoutHeap_DFS(P, H_G, 0, G)
    
    return empty_seq


class EHeapElement:
    def __init__(self, seq, weight):
        self.seq = seq
        self.weight = weight
    def __lt__(self, other):
        return self.weight < other.weight
    # def __str__(self):
    #     return f"{self.seq}"
    def __repr__(self):
        return f"EHE({self.weight}'\t'{self.seq})"
    def __hash__(self):
        return hash((tuple(self.seq), self.weight))
    def __eq__(self, other) : 
        return self.__dict__ == other.__dict__
    def __ne__(self, other):
        if not(isinstance(other, EHeapElement)): return True
        if not(self.__class__ == other.__class__): return True
        return not(self.__dict__ == other.__dict__)


# node: roote node to start BFS from
def P_to_Heap(H, P, node):
    Er = EHeapElement([node], 0)

    # if the empty sequence is the only sequence for a path from s to t
    if P.adj[node] == {}:
        H.add_node(Er)
        return Er

    n1 = next(iter(P.adj[node]))
    En = EHeapElement([node, n1], P.adj[node][n1]['weight'])
    H.add_edge(Er, En)

    # BFS P
    nx.set_node_attributes(P, False, 'visited')
    queue = []
    w_t = P.adj[node][n1]['weight'] # total weight
    queue.append((n1, Er, w_t))

    while queue:
        l = queue.pop(0)
        l_t = l[0]      # tail
        l_s = l[1].seq
        w_t = l[2]

        # new sequence (already added but we need to construct the parent)
        new_seq = l_s + [l_t]
        p = EHeapElement(new_seq, w_t)

        for child in P.adj[l_t]:
            if not P.nodes[child]['visited']:
                P.nodes[child]['visited'] = True
                if P.adj[l_t][child]['cross_edge'] == True:
                    # handle child
                    t_w = w_t + P.adj[l_t][child]['weight'] # total weight
                    c = EHeapElement(new_seq + [child], t_w)

                    # add parent-child connection to the heap
                    H.add_edge(p, c)

                    # add child to queue
                    queue.append( (child, p, t_w) )
                else:
                    t_w = w_t + P.adj[l_t][child]['weight'] # total weight
                    c = EHeapElement(l_s + [child], t_w )   # add child to sequence
                    H.add_edge(p, c)
                    queue.append( (child, l[1], t_w) )
    return Er

# root: EHeapElement consisting of a sequence of STCedges and the additional stc
# H: The heap to pop from
def pop_from_H(root, H):

    # check if there is anything to pop
    if H.adj[root]:
        new_root = min(H.adj[root])
    else:
        return root, None
    
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
    return v, new_root


def get_paths_between(head, tail, pred, path=[]):
    path = path + [head]
    if head == tail:
        return [path]
    if pred[head] == []:
        return []
    paths = []
    for node in pred[head]:
        newpaths = get_paths_between(node, tail, pred, path)
        for newpath in newpaths:
            paths.append(newpath)
    return paths

# for each pair of subsequent STCedges in the sequence, we want to find all
# paths that connect them, consisting of non-STCedges only.
# For each of these paths, we want to continue doing this in the following pairs
# of subsequent STCedges
# p, EHeapElement: contains sequence of sidetrack edges and additional cost
# sd: shortest path distance, to be added to the additional cost of a path
# pred: predecessor nodes on the shortest path tree
# dst: the destination node
def get_all_paths_for_sequence(p, sd, pred, dst):
    all_paths = [[]] # all paths corresponding to p
    seq = p.seq
    seq.append(STCedge(dst, None, {'weight': 0, 'sidetrackCost': 0}))

    for i in range(len(seq)):
        # find all shortest paths from e.head to next(e).tail
        if i < len(seq) - 1:
            new_all_paths = []
            m = get_paths_between(seq[i].head, seq[i+1].tail, pred)
            for l in all_paths:
                for n in m:
                    new_all_paths.append(l + n)
            all_paths = new_all_paths

    # cost of the entire path
    c = sd + p.weight
    all_paths = [(sb, c) for sb in all_paths]
    return all_paths


# G:    a networkx DiGraph
# src:  the source node
# dst:  the destination node
# k:    the number of paths
def k_shortest_paths(G, src, dst, k):
    st_all = time.time()
    # Calculate single-destination shortest path tree; this is the
    # same as a single-source shortest path tree in G-reverse
    st_rev = time.time()
    R = G.reverse(copy=True)
    et_rev = time.time()

    # Essentially calculate shortest path tree but just remember
    # predecessor of node on shortest path to dst and min distance to each node from
    # a given source
    st_dijk = time.time()
    pred, dist = nx.dijkstra_predecessor_and_distance(R, dst) # ( should take O(m + nlogn) )
    et_dijk = time.time()

    # Check if src and dst are connected
    if dist.get(src) == None:
        return []

    # remove nodes from which dst can't be reached
    clean_graph(G, dist)

    # Calculate sidetrack costs for every edge and add them as attribute to the edge
    st_calc_sidetrack_cost = time.time()
    calc_sidetrack_cost(G, dist)
    et_calc_sidetrack_cost = time.time()

    # Calculate Hout(v) for every vertex and add the resulting heap as attribute to v
    st_calc_Hout = time.time()
    for v in G.nodes():
        G.nodes[v]['Hout'] = Hout(G,v)
    et_calc_Hout = time.time()

    # Calculate H_G(v) for every vertex, a balanced heap of the roots of Hout on the path from v to t
    # where every root also points to the rest of its original Hout heap
    st_calc_HG = time.time()
    G = calc_H_G(G, pred, dst) # For every vertex H_G is added as an attribute
    et_calc_HG = time.time()

    P = nx.DiGraph()
    st_prepare_and_augmentP = time.time()
    Proot = prepare_and_augmentP(P, G, src) # results in a graph of STCedge elements
    et_prepare_and_augmentP = time.time()
    # P is now the completed path graph P(G), its root is returned by the function above

    # Lastly, P(G) will need to be transformed into a 4-heap H(G), so that the nodes
    # in H(G) represent paths in G. H(G) is constructed by forming a node for each path in
    # P(G) rooted at Proot.
    H = nx.DiGraph()
    st_P_to_Heap = time.time()
    Hroot = P_to_Heap(H, P, Proot)
    et_P_to_Heap = time.time()

    # Convert STCedge sequence to sequence of nodes describing a path in G
    # Convert cost to complete cost
    # Goal: a list of tuples (l, c): l, a list of nodes; c, the cost
    paths = []
    sd = dist[src] # shortest path distance, to be added to the sidetrack cost

    # To find the k shortest paths, pop at most k times from H and append to list
    # We will convert STCedge sequences to a list of paths after each pop
    last = 0
    st_popping = time.time()
    while len(paths) < k:

        if Hroot is None:
            if len(paths) < k:
                print(f"No more than {len(paths)} paths")
            break

        p, Hroot = pop_from_H(Hroot, H)
        paths += get_all_paths_for_sequence(p, sd, pred, dst)

        if Hroot is None:
            if len(paths) < k:
                print(f"No more than {len(paths)} paths")
            break

        while True and Hroot.weight == last: # pop all sequences with same weight, for testing purposes, can be set to False for better performance

            last = Hroot.weight
            p, Hroot = pop_from_H(Hroot, H)
            paths += get_all_paths_for_sequence(p, sd, pred, dst)

            if Hroot is None:
                if len(paths) < k:
                    print(f"No more than {len(paths)} paths")
                break
        if Hroot != None:
            last = Hroot.weight

    et_popping = time.time()
    et_all = time.time()

    t_rev = et_rev - st_rev
    t_dijk = et_dijk - st_dijk
    t_calc_sidetrack_cost = et_calc_sidetrack_cost - st_calc_sidetrack_cost
    t_calc_Hout = et_calc_Hout - st_calc_Hout
    t_calc_HG = et_calc_HG - st_calc_HG
    t_prepare_and_augmentP = et_prepare_and_augmentP - st_prepare_and_augmentP
    t_P_to_Heap = et_P_to_Heap - st_P_to_Heap
    t_popping = et_popping - st_popping
    t_all = et_all - st_all

    t_l = [('t_rev', t_rev),
           ('t_dijk', t_dijk),
           ('t_calc_sidetrack_cost', t_calc_sidetrack_cost),
           ('t_calc_Hout', t_calc_Hout),
           ('t_calc_HG', t_calc_HG),
           ('t_prepare_and_augmentP', t_prepare_and_augmentP),
           ('t_P_to_Heap', t_P_to_Heap),
           ('t_popping', t_popping)
           ]
    for name, t in t_l:
        print(f"{round(t/t_all * 100)}% {name}")


    return paths


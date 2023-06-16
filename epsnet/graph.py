# Code from eps-net repository by Patrick Emonts 2023

import networkx as nx
import numpy as np
from heapq import heappush, heappop
import copy
import sys
from collections import defaultdict
from itertools import islice
import logging


def generate_random_graph(n,p,weight_interval=[0,10]):
    """Generate a random graph

    Args:
        n (int): Number of edges
        p (float): Probability of an edge
        weight_interval (tuple, optional): Range of weights. Defaults to [0,10].

    Returns:
        nx.DiGraph: Random digraph with weighted edges
    """
    # Generate the graph ourselves
    g_rand = nx.DiGraph()
    g_rand.add_nodes_from(range(n))

    for u in g_rand.nodes:
        for v in range(u+1, len(g_rand.nodes)):
            # With prob 0.5, we add an edge
            if np.random.uniform() < p:
                weight = np.random.randint(*weight_interval)
                g_rand.add_edge(u, v, weight=weight)
    return g_rand


def next_best_paths_waterman(graph, source, target, cutoff,print_progress=False):
    """Compute the next best path between source and target in graph.
    The cut-off is given as an absolute weight of the path.
    The algorithm by Waterman and Byres (10.1016/0025-5564(85)90096-3) is used to find the next best.

    Args:
        graph (nx.DiGraph): Directed Graph
        source (label): Label of the source node
        target (label): Label of the target node
        cutoff (float): Weight used as a cut-off for the next-best paths
        print_progress (bool, optional): Print the progress during the algorithm. Defaults to False.

    Returns:
        list: list of paths through the graph
    """
    shortest_path_dict = nx.shortest_path_length(graph,target=target,weight="weight")
    # If a given node has no outgoing paths, we set the cost to infinity
    for node in graph.nodes:
        if node not in shortest_path_dict:
            shortest_path_dict[node]=np.inf
    # We use a list as a stack
    # We push 3-tuples with the information (current node, next node, distance of next node from source)
    stack = []
    dest = []
    # Fill the first entry of the stack
    source_successors = graph.successors(source)
    for n in source_successors:
        weight = graph.get_edge_data(source,n)["weight"]
        if weight + shortest_path_dict[n]<cutoff:
            stack.append((source,n,weight,[source]))

    def rec_search(graph, source, target, cutoff, stack):
        # The path tracking is still problematic here!
        if len(stack)>0:
            node, next_node, d, path= stack.pop()
            if print_progress:
                print(node, next_node, d, path)
            path_loc = path.copy()
            path_loc.append(next_node)
            if next_node == target:
                dest.append((path_loc,d))
            else:
                successors = graph.successors(next_node)
                for n in successors:
                    weight = graph.get_edge_data(next_node,n)["weight"]
                    if weight + shortest_path_dict[n]+d < cutoff:
                        stack.append((next_node,n,weight+d,path_loc))
            rec_search(graph, source, target, cutoff, stack.copy())

    rec_search(graph,source,target,cutoff, stack)

    return sorted(dest,key=lambda x: (x[1],tuple(x[0]))) 

def next_best_k_paths_waterman(graph, source, target, k ,print_progress=False):
    """Compute the next best path between source and target in graph.
    The cut-off is given in percent of the total cost.
    The output is guranateed to be of length k. If multiple paths of the same length exist, the cut might be arbitrary.

    The algorithm by Waterman and Byres (10.1016/0025-5564(85)90096-3) is used to find the next best.

    Args:
        graph (nx.DiGraph): Directed Graph
        source (label): Label of the source node
        target (label): Label of the target node
        k (int): number of next-best paths
        print_progress (bool, optional): Print the progress during the algorithm. Defaults to False.

    Returns:
        list: list of paths through the graph
    """
    min_path = nx.shortest_path(graph,source,target,weight="weight")
    nsol = 1
    solutionvec = [min_path]
    min_weight = nx.path_weight(graph,min_path,"weight")
    cutoff = min_weight
    while(nsol<k):
        increment = max(0.2, 1.1*cutoff)
        cutoff += increment
        next_best_solvec=next_best_paths_waterman(graph, source, target, cutoff)
        nsol = len(next_best_solvec)
    return next_best_solvec[:k]


def next_best_k_paths_dijkstra(graph,source, target, k):
    path_set = set()
    countdict = {}
    for u in graph.nodes():
        countdict[u]=0
    path_heap = [] # to be used as heap
    heappush(path_heap,(0,[source])) #insert path Ps = {s} into B with cost 0
    while len(path_heap)>0 and countdict[target]<k: #while B is not empty and countt < K:
        cost_u, path_u = heappop(path_heap) #– let Pu be the shortest cost path in B with cost C
        u = path_u[-1]
        countdict[u]+=1 #– B = B − {Pu }, countu = countu + 1
        if u == target:
            path_set.add((cost_u,tuple(path_u)))
            #– if u = t then P = P U {Pu}
        if countdict[u]<k: #– if countu ≤ K then
            for v in graph.neighbors(u): #for each vertex v adjacent to u:
                path_v = path_u +[v]
                cost_v = cost_u+graph.edges[u,v]["weight"]# let Pv be a new path with cost C + w(u, v) formed by concatenating edge (u, v) to path Pu
                heappush(path_heap,(cost_v,path_v)) #– insert Pv into B
    # Invert the tuples to match convention
    dest = [(list(path),cost) for cost,path in sorted(list(path_set))]
    return dest


class DataTuple:
    def __init__(self,id,data):
        self.id = id
        self.data = data

    def __lt__(self, other):
        return self.id<other.id

class EppsteinShortestPathAlgorithm:
    
    """ Adapted from https://github.com/lppcom/Eppstein-Algorithm-in-Python"""

    def __init__(self, graph, source ='s', destination ='t'):
        self._G = graph
        self.path_tree = []
        self.source = source
        self.destination = destination
        self.sidetrack_edges=[]
        self.shortest_path_distance = None
        self.counter = 0
        self._pre_process()

    def _get_path_from_predecessors(self,pred={}, destination=None):
        if not destination:
            destination = self.destination
        path = list()
        while(True):
            value = pred.get(destination, None)
            if not value:
                break
            else:
                # Get rid of the list
                value = value[0]
            path.append((value, destination))
            destination = value
        path = list(reversed(path)) # reverse the path
        return path

    def _all_shortest_paths(self):
    
        """ Find all shortest paths from every node to destination """
        #make a reversed graph (reversing all the edges), so we can find single destination shortest paths problem
        
        _reverse_graph =  self._G.reverse(copy=True)
        _reverse_pred, _dist = nx.bellman_ford_predecessor_and_distance(_reverse_graph,self.destination) 
        logging.debug("EppsteinShortestPath: reverse_pred, & dist by using bellman ford ")
        _pred = defaultdict(dict)
        for node, neighbor in _reverse_pred.items():
            # We assume here that we have one predecessor
            if neighbor is not None and len(neighbor) > 0:
                _pred[neighbor[0]]=node
        for counter, node in enumerate(self._G.nodes()):
            try:
                self._G.nodes[node]['target_distance']=_dist[node]
            except KeyError:
                self._G.nodes[node]['target_distance']=float('inf')
            _path=self._get_path_from_predecessors((_reverse_pred), destination=node)
            path = list(reversed([(value, key) for key,value in _path]))
            self._G.nodes[node]['path']=path
            
        self.shortest_path_distance = self._G.nodes[self.source]['target_distance']


    def __get_sidetrack_edges(self, G):
        sp_edges = set([i for i in (path for node in G.nodes() for path in G.nodes[node]['path'])])
        all_edges = set([edge for  edge in G.edges()])
        sidetrack_edges = all_edges - sp_edges
        return sidetrack_edges

    def _init_path_tree(self):
        sp = self._G.nodes[self.source]['path']
        sp_distance = self._G.nodes[self.source]['target_distance']
        node_info = {}
        node_info['sigma_e']=0
        node_info['path']=sp
        node_info['distance']=sp_distance
        node_info['visited_edges']=set()
        node_info['source']=None
                
        heappush(self.path_tree, DataTuple(node_info['sigma_e'], node_info))
        
    def _head_tail(self, edge):
        return edge[1], edge[0]
    

    def _is_tree(self, G):
        if nx.number_of_nodes(G) != nx.number_of_edges(G)+1:
            return False
        return nx.is_connected(G)
    
    def _build_graph(self, adj_dict={}, is_directed = True):
        """
        Build Graph from supplied adjacency list or adj dict
        """
        if is_directed:
            G = nx.DiGraph()
        else:
            G = nx.Graph()
        for node in adj_dict.keys():
            G.add_node(node)
        for node, neighbor_list in adj_dict.items():
            for neighbor in neighbor_list:
                G.add_edge(node, neighbor)
        return G
        
    def _has_path(self, source=None, destination=None, edges=[]):
        """
        from the list of given edges determine whether a path exists in betweeen source and destination
        """
        if source is None:
            source = self.source
        if destination is None:
            destination = self.destination
        adj_dict = self._adj_dict(edges)

        G = self._build_graph(adj_dict)
        return nx.has_path(G, source, destination)
    
    def _adj_dict(self, edges):
        """
        returns adjacency dict for given set of edges
        """
        adj_dict = defaultdict(list)
        for i, j in edges:
            adj_dict[i].append(j)
        return adj_dict
    
    def __is_valid_sidetrack_edge(self, path, edge, ):
        # Check whether the given edge can become a side track edge for given path
        if not path:
            return False
        if edge in path:
            # logging.debug('Edge is in path' + str(path) +", Edge:"+str(edge))
            return False
        adj_dict = self._adj_dict(path)
        head, tail = self._head_tail(edge)
        
        if adj_dict[tail]:
            return True
        return False
    
    def __get_sidetrack_path(self, path, sidetrack_edge):
        head, tail = self._head_tail(sidetrack_edge)
        to_remove = []
        for counter, pe in enumerate(path):
            ph,pt = self._head_tail(pe)
            if pt == tail:
                to_remove.append(counter)
        # make sure that we are dealing with a path
        assert (len(to_remove) == 1)
        remaining_edges = [pe for counter, pe in enumerate(path) if counter not in to_remove]
        remaining_edges.insert(to_remove[0], sidetrack_edge)
        remaining_edges.extend(list(self._G.nodes[head]['path']))
        return self._get_path_from_edges(remaining_edges)
    
    def _get_path_from_edges(self, edges, destination=None):
        if destination is None:
            destination = self.destination
        if self._has_path(edges=edges, source=self.source, destination = destination):
            adj_dict = self._adj_dict(edges)
            G = self._build_graph(adj_dict = adj_dict)
            pred, dist = nx.bellman_ford_predecessor_and_distance(G, self.source, destination)
            return self._get_path_from_predecessors(pred, destination)
    
    def __get_sidetrack_edge_info(self, edge=None, prev_sigma_e=None, path=None, source = None):
        info = {}
        head,tail = self._head_tail(edge)
        info['sigma_e'] = self._G.nodes[head]['target_distance']-self._G.nodes[tail]['target_distance']+self._G.edges[tail,head]['weight']+prev_sigma_e
        info['edge'] = edge
        info['source'] = source
        info['path'] = self.__get_sidetrack_path(path, edge)
        return info
    
    def _insert_edge(self, root, children=None, node_info=None):
        """
        draw an edge in between root and the given node
        """
        graph = self.path_tree
        graph.add_node(children, index=children, node_info=node_info)
        graph.add_edge(root, children)
        

    def _build_path_tree(self,source=None, path=[], sidetrack_edges=set(), prev_sigma_e=0, current_vertex=0,visited_edges=set()):
        if source is None:
            source = self.source
        for edge in sidetrack_edges:
            if self.__is_valid_sidetrack_edge(path=path, edge=edge, ):
                info = self.__get_sidetrack_edge_info(edge=edge, source=source, prev_sigma_e=prev_sigma_e, path=path)
                node_info = {}
                node_info['prev_sigma_e'] = prev_sigma_e
                node_info['sigma_e'] = info['sigma_e']
                node_info['path'] = info['path']
                node_info['distance'] = self.shortest_path_distance + node_info['sigma_e']
                seen_edges = copy.deepcopy(visited_edges)
                seen_edges.add(edge)
                node_info['visited_edges'] = visited_edges
                node_info['source'] = info['edge']
                flag = True
                for ancestor_edge in visited_edges:
                    if ancestor_edge not in info['path']:
                        flag = False
                if flag:
                    heappush(self.path_tree, DataTuple(node_info['sigma_e'], node_info))
            
    
    def _pre_process(self):
        """
            Finall all shortest path from each node to destination.
            Retrieve side track edges
            build side_track
        """
        self._all_shortest_paths()
        logging.info("EppsteinShortestPath: all shortest paths completed")
        sidetrack_edges = self.__get_sidetrack_edges(self._G)
        self.sidetrack_edges = sidetrack_edges
        self._init_path_tree()
        logging.info("EppsteinShortestPath: Initialization has been done")
        
    def retrieve_k_best(self):
        while(self.path_tree):
            data_el_info = heappop(self.path_tree)
            el_info = data_el_info.data
            visited_edges = el_info['visited_edges']
            seen_edges = copy.deepcopy(visited_edges)
            seen_edges.add(el_info['source'])
            
            sigma_e = el_info['sigma_e']
            distance = el_info['distance']
            path = el_info['path']
            yield path, distance
            source = el_info['source']
            if data_el_info.id==0:
                self._build_path_tree(sidetrack_edges = self.sidetrack_edges, path=path)
            else:
                self._build_path_tree(source = source,path = path, sidetrack_edges = self.sidetrack_edges, prev_sigma_e = sigma_e,visited_edges=seen_edges)

    def get_successive_shortest_paths(self, vertex_mode=True):
        for path, weight in self.retrieve_k_best():
            if vertex_mode:
                # Convert edges in the convention (A,B),(B,C) to a list of vertices [A,B,C]
                vertex_path = []
                for ind, edge in enumerate(path):
                    if ind == 0:
                        vertex_path.append(edge[0])
                    vertex_path.append(edge[1])
                yield vertex_path, weight
            else:
                yield path, weight
        
    def get_k_shortest_paths(self,k:int):
        """Return the k shortest paths, sorted first by their weight and then by their path considered as a tuple

        Args:
            k (int): Number of paths

        Returns:
            list: List of (path,weight) tuples
        """
        pathvec = list(islice(self.get_successive_shortest_paths(),k))
        return sorted(pathvec, key=lambda x:(x[1],tuple(x[0])))
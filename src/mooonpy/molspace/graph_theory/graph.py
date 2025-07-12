# -*- coding: utf-8 -*-
"""
This module provides basic support generating and working with graphs.
The main purpose is to provide in interface between a Molspace instance
and the graph data structure used in all the graph theory supported 
workflows.
"""
from ..molspace import Molspace


def generate_graph(mol: Molspace):
    """
    Generates an undirected graph from a Molspace instance.
    
    .. note::
        The generated graph will contain the most current nodes/edges each
        time this function is called, so if adding atoms/bonds to a Molspace
        instance, it is important to first update the Molspace instance before
        calling this function.

    :param mol: Molspace instance to generate undirected graph from
    :type mol: Molspace
    :return: graph
    :rtype: dict[int, list[int]]
    """
    graph = {i:[] for i in mol.atoms} # {atomID : [bonded-neighbors] }
    for (id1, id2) in mol.bonds:
        graph[id1].append(id2)
        graph[id2].append(id1)
    return graph

            
            
def find_cumulative_neighs(graph: dict[int, list[int]], node: int, max_depth: int):
    """
    Finds the cummulative neighbors of a graph starting at
    a given node and traversing to a maximum depth of max_depth.
    Depth is defined as follows:
        
    1. first neighbors
    2. second neighbors
    3. thrid neighbors
    
    .. note::
        For rings/cycles, where each neighbor has two classifications of
        depth due to ring symmetry, the lowest neighbor depth classification
        will be assigned and it will not be duplicated in the higher neighbor
        depth set.
        

    :param graph: An undirected graph
    :type graph: dict[int, list[int]]
    :param node: A node to find all cummulative neighbors from
    :type node: int
    :param max_depth: The maximum traversal depth to search (typically 4)
    :type max_depth: int
    :return: neighbors
    :rtype: dict[int, set[int]]
    """
    neighbors = {i+1 : set() for i in range(max_depth)} # { depth : {id1, id2, ...} }
    neighbors[1] = set(graph[node]) 
    visited = set(graph[node] + [node])
    for depth in range(1, max_depth):
        for i in neighbors[depth]:
            for j in graph[i]:
                if j in visited: continue  
                neighbors[depth+1].add(j)
                visited.add(j)
    return neighbors
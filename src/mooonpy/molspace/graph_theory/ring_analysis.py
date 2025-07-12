# -*- coding: utf-8 -*-
"""
This module provides basic support for finding cycles/rings in
a molecular system.
"""
from .graph import find_cumulative_neighs
import math


def _reduced_graph(graph, node, max_depth):
    # We need one additionaly neighbor to generate a graph that terminates after
    # the max_depth cutoff, hence the "max_depth+1" in find_cumulative_neighs()    
    neighs = find_cumulative_neighs(graph, node, max_depth+1)    
    rgraph = {node : graph[node]} # {atomID : [bonded-neighbors] }
    for depth in neighs: 
        for i in neighs[depth]:
            if depth < len(neighs): 
                rgraph[i] = graph[i]
            else: rgraph[i] = [] # graph termination
    return rgraph, neighs


def _dfs_cycles(graph, start, max_ring=None, visited=None):
    if visited is None: visited = set()
    if max_ring is None: max_ring = len(graph)
    stack = [(start, [])]
    while stack:
        node, path = stack.pop()
        if node == start and len(path) > 2: 
            # Canonical form (e.g., [1,2,3] == [2,3,1])
            if path[0] > path[-1]:
                path = path[::-1]
            min_index = path.index(min(path))
            canonical = tuple(path[min_index:] + path[:min_index])
            visited.add(node)
            yield canonical
        if node in visited: continue
        for neighbor in graph[node]:
            if neighbor in path: continue
            new_path = path + [neighbor]
            if len(new_path) <= max_ring:
                stack.append((neighbor, new_path))
            

def find_rings(graph: dict[int, list[int]], ring_sizes: tuple[int]=(3,4,5,6,7)):
    """
    Finds all cycles/rings in a graph that are specified in the ring_sizes
    parameter.
    
    .. note::
        Each ring will be sorted in the order it was traversed in the graph and
        will be in the canonical form (e.g., canonical=[1,2,3] vs rotated=[2,3,1]).

    :param graph: An undirected graph
    :type graph: dict[int, list[int]]
    :param ring_sizes: Tuple containing ring sizes to search for in graph
    :type ring_sizes: tuple[int]
    :return: rings
    :rtype: list[tuple[int]]
    """
    # Find the maximum depth we need to search based on the largest ring we expect. 
    # We only need to search a depth of half of the largest ring we want to find 
    # since the ring is symmetric about that many neighbors. It is also quicker to
    # use the "in" keyword on a set, so convert ring_sizes to a set for performance
    max_ring = max(ring_sizes)
    max_depth = math.ceil(0.5*(max_ring))
    rings2check = set(ring_sizes)


    # Start walking around the graph to find the rings. We will need to keep
    # track of rings and permuations of the atomIDs in each ring that where
    # already walked along, which will be accomplished by sorting the atomIDs
    # and logging the ring into sorted_rings. The return value will be unique
    # rings where the atomIDs are ordered in the direction they where walked.
    sorted_rings = set()
    visited = set()
    rings = set()
    for node in graph:  
        if len(graph[node]) <= 1 or node in visited: continue

        # Reduce graph for performance to smallest possible
        # adjacency list for quickest possible dfs traversal
        rgraph, neighs = _reduced_graph(graph, node, max_depth)
        
        # Use a depth first search traveral on the reduced graph
        for path in _dfs_cycles(rgraph, node, max_ring=max_ring, visited=visited):
            sorted_ring = tuple(sorted(path))
            if len(sorted_ring) in rings2check and sorted_ring not in sorted_rings:                
                rings.add( path )
                sorted_rings.add( sorted_ring )
        
        # Add node to visted so it is not walked along again in the next dfs
        visited.add(node)
        
    return sorted(rings, key=lambda x: min(x))
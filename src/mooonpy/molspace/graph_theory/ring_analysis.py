# -*- coding: utf-8 -*-
"""
This module provides basic support for finding cycles/rings in
an undirected graph.
"""
from .interface import find_cumulative_neighs
import math


def _reduced_graph(graph, node, max_depth):
    # We need one additionaly neighbor to generate a graph that terminates after
    # the max_depth cutoff, hence the "max_depth+1" in find_cumulative_neighs()    
    neighs = find_cumulative_neighs(graph, node, max_depth + 1)
    rgraph = {node: graph[node]}  # {atomID : [bonded-neighbors] }
    for depth in neighs:
        for i in neighs[depth]:
            if depth < len(neighs):
                rgraph[i] = graph[i]
            else:
                rgraph[i] = []  # graph termination
    return rgraph


def _dfs_cycles(graph, start, max_ring):
    stack = [(start, [start], {start})]
    while stack:
        current, path, visited = stack.pop()
        path_size = len(path)
        for neighbor in graph[current]:
            if neighbor == start and path_size > 2:
                yield path
            elif neighbor not in visited and neighbor >= start:
                new_path = path + [neighbor]
                if len(new_path) > max_ring: continue
                stack.append((neighbor, new_path, visited | {neighbor}))


def ring_canonical(ring: list) -> tuple:
    """
    Sort list of ring ID's to canonical form of "lowest ID, lower neighbor-wise"
    ie [5,254,6,12,4,3] -> (3,4,12,6,254,5)
    with 3 as lowest, and 4 as the lower of 4 and 5, reversing the direction
    """
    min_index = ring.index(min(ring))
    if ring[min_index + 1] < ring[min_index - 1]:
        canonical = tuple(ring[min_index:] + ring[:min_index])  # forward shift
    else:
        ring.reverse()
        min_index = len(ring) - min_index - 1
        canonical = tuple(ring[min_index:] + ring[:min_index])  # reverse then shift
    return canonical


def find_rings(graph: dict[int, list[int]], ring_sizes: tuple[int] = (3, 4, 5, 6, 7)):
    """
    Finds all cycles/rings in a graph that are specified in the ring_sizes
    parameter.
    
    .. note::
        Each ring will be sorted in the order it was traversed in the graph and
        will be in the canonical form (e.g., canonical=(1,2,3) vs non_canonical=(2,3,1)).

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
    max_depth = math.ceil(0.5 * (max_ring))
    rings2check = set(ring_sizes)

    # Start walking around the graph to find the rings. We will need to keep
    # track of rings and permutations of the atomIDs in each ring that where
    # already walked along, which will be accomplished by sorting the atomIDs
    # and logging the ring into sorted_rings. The return value will be unique
    # rings where the atomIDs are ordered in the direction they were walked.
    sorted_rings = set()
    rings = set()
    for node in graph:
        if len(graph[node]) <= 1: continue

        # Reduce graph for performance to smallest possible
        # adjacency list for quickest possible dfs traversal
        rgraph = _reduced_graph(graph, node, max_depth)

        # Use a depth first search traversal on the reduced graph
        for ring in _dfs_cycles(rgraph, node, max_ring):
            sorted_ring = tuple(sorted(ring))
            ring_size = len(ring)
            if ring_size in rings2check and sorted_ring not in sorted_rings:
                canonical = ring_canonical(ring)

                # print(ring, canonical)
                rings.add(canonical)
                sorted_rings.add(sorted_ring)

    return sorted(rings, key=lambda x: min(x))

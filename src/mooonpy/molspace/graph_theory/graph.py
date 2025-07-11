# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 13:10:48 2025

@author: jdkem
"""
import math


def generate_graph(mol):
    graph = {i:[] for i in mol.atoms} # {atomID : [bonded-neighbors] }
    for (id1, id2) in mol.bonds:
        graph[id1].append(id2)
        graph[id2].append(id1)
    return graph



            
            
def find_cumulative_neighs(graph, node, max_depth):
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

    

def reduced_graph(graph, node, max_depth):
    # We need one additionaly neighbor to generate a graph that terminates after
    # the max_depth cutoff, hence the "max_depth+1" in find_cumulative_neighs()    
    neighs = find_cumulative_neighs(graph, node, max_depth+1)    
    rgraph = {node : graph[node]} # {atomID : [bonded-neighbors] }
    for depth in neighs: 
        for i in neighs[depth]:
            if depth < len(neighs): 
                rgraph[i] = graph[i]
            else: rgraph[i] = [] # graph termination
    return rgraph


def dfs_cycles(graph, start, end, visited=None):
    if visited is None: visited = set()
    stack = [(start, [])]
    while stack:
        node, path = stack.pop()
        if len(path) >= 3 and node == end: 
            visited.add(node)
            yield path
            continue
        if node in visited: continue
        for neighbor in graph[node]:
            if neighbor in path: continue
            stack.append((neighbor, path+[neighbor]))
            

def find_rings(graph, ring_sizes=(3,4,5,6,7)):
    # Find the maximum depth we need to search based on the largest ring we 
    # expect. We only need to search a depth of half of the largest ring we
    # want to find as the ring is symmetric about that many neighbors. It is
    # also quicker the use the "in" keyword on a set, so convert ringsizes
    # to a set for quicker speed
    max_depth = math.ceil(0.5*(max(ring_sizes)))
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
        if len(graph[node]) <= 1: continue
        
        # Reduce graph for performance to smallest possible
        # adjacency list for quickest possible dfs traversal
        rgraph = reduced_graph(graph, node, max_depth)
        
        # Use a depth first search traveral on the reduced
        # graph and find all cycles in reduced graph
        for path in dfs_cycles(rgraph, start=node, end=node, visited=visited):
            sorted_ring = tuple(sorted(path))
            if len(sorted_ring) in rings2check and sorted_ring not in sorted_rings:
                if path[0] < path[-1]:
                    ordered_ring = tuple(path)
                else: ordered_ring = tuple(path[::-1])
                rings.add( ordered_ring )
                sorted_rings.add( sorted_ring )
        visited.add(node)
    return rings
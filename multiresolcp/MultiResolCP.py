import _km_omp as _cp
import networkx as nx
from itertools import compress
import numpy as np
import scipy
from scipy.sparse import triu 

def detect(G, nodes_in_part1, nodes_in_part2, part_to_project, resol = 1, node_capacity = {}, num_samples = 100, num_runs=10, consensus_threshold=0.9, significance_level = 0.05, num_rand_nets = 500):

	# Input check 
	if part_to_project == 'part1':
		nodes_to_project = nodes_in_part1
		nodes_to_be_collapsed = nodes_in_part2
		pass	
	elif part_to_project == 'part2':
		nodes_to_project = nodes_in_part2
		nodes_to_be_collapsed = nodes_in_part1
	else:
		raise Exception("Error in part_to_project. Set either part_to_project = 'part1' or part_to_project = 'part2.'") 

	nodes_to_project = set(nodes_to_project)	
	nodes_to_be_collapsed = set(nodes_to_be_collapsed)	
	if len(list(nodes_to_project.intersection(nodes_to_be_collapsed))) >0:
		raise Exception("Error in nodes_in_part1 and nodes_in_part2. nodes_in_part1 and nodes_in_part2 should not contain the same node.") 

	if nodes_to_project.union(nodes_to_be_collapsed)  != set(G.nodes()):
		raise Exception("Error in nodes_in_part1 and nodes_in_part2. Some nodes are missing in nodes_in_part1 and nodes_in_part2.") 
		raise Exception("Error in part1 and nodes_to_be_collapsed. ") 
		
	if len(node_capacity) == 0:
		node_capacity = np.array(np.ones(len(nodes_to_be_collapsed))).astype(float)
		nodes_to_project = list(nodes_to_project)
		nodes_to_be_collapsed = list(nodes_to_be_collapsed)
	else:
		if len(set(node_capacity.keys()).symmetric_difference(nodes_to_be_collapsed))>9:
			raise Exception("Error in node_capacity. Some nodes are missing in node_capacity.") 
		nodes_to_project = list(nodes_to_project)
		nodes_to_be_collapsed = list(nodes_to_be_collapsed)
		node_capacity = np.array([node_capacity[r] for r in nodes_to_be_collapsed]).astype(float)

	Np = len(nodes_to_project)
	Nr = len(nodes_to_be_collapsed)

	# Make the list of edges in the given network
	A = nx.adjacency_matrix(G, nodes_to_project + nodes_to_be_collapsed)
	r, c = triu(A).nonzero()
	edges = np.array([[rc[0], rc[1]] for rc in zip(r, c)]).astype(int)
	
	# Pass parameters to a c++ function (src/_km_ompnet.cpp)	
	results = _cp._detect(edges = edges,\
			nodes_to_project = np.array(list(range(Np))).astype(int),\
			nodes_to_be_collapsed = np.array(list(range(Np, Nr + Np))).astype(int),\
			node_capacity = np.array(node_capacity).astype(float),\
			resol = float(resol),\
			num_samples = int(num_samples),\
			num_runs = int(num_runs),\
			consensus_threshold = float(consensus_threshold),\
			significance_level = float(significance_level),\
			num_rand_nets = int(num_rand_nets))

	# Retrieve the results 
	c = results[0].astype(int)	
	x = results[1]
	
	# Remove the homeless nodes that do not belong to any consensus CP pairs  	
	b = c>=0
	c = dict(zip(compress(nodes_to_project, b), c[b]))	
	x = dict(zip(compress(nodes_to_project, b), x[b]))
	
	return c, x	

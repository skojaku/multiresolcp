import _km_omp as _cp
import networkx as nx
from itertools import compress
import numpy as np
import scipy
from scipy.sparse import triu 

def detect(G, nodes_in_part1, nodes_in_part2, part_to_project, resol = 1, node_capacity = {}, num_samples = 100, consensus_threshold=0.9, significance_level = 0.05, num_rand_nets = 500):


	""" INPUT VALIDATION """ 
	if part_to_project == 'part1':
		_nodes_side_A = nodes_in_part1
		_nodes_side_B = nodes_in_part2
		pass	
	elif part_to_project == 'part2':
		_nodes_side_A = nodes_in_part2
		_nodes_side_B = nodes_in_part1
	else:
		raise Exception("Invalid input part_to_project. Set either part_to_project = 'part1' or part_to_project = 'part2'.") 

	_nodes_side_A = set(_nodes_side_A)	
	_nodes_side_B = set(_nodes_side_B)	
	if len(list(_nodes_side_A.intersection(_nodes_side_B))) >0:
		raise Exception("Invalid inputs nodes_in_part1 and nodes_in_part2. nodes_in_part1 and nodes_in_part2 should not contain the same node.") 

	if _nodes_side_A.union(_nodes_side_B)  != set(G.nodes()):
		raise Exception("Invalid inputs nodes_in_part1 and nodes_in_part2. Some nodes are missing.") 
		
	if len(node_capacity) == 0:
		node_capacity = np.array(np.ones(len(_nodes_side_B))).astype(float)
		_nodes_side_A = list(_nodes_side_A)
		_nodes_side_B = list(_nodes_side_B)
	else:
		if len(set(node_capacity.keys()).symmetric_difference(_nodes_side_B))>0:
			raise Exception("Invalid input node_capacity. Some nodes are missing in node_capacity.") 
		_nodes_side_A = list(_nodes_side_A)
		_nodes_side_B = list(_nodes_side_B)
		node_capacity = np.array([node_capacity[r] for r in _nodes_side_B]).astype(float)


	""" CORE-PERIPHERY DETECTION  """ 
	# Make the list of edges in the given network
	A = nx.adjacency_matrix(G, _nodes_side_A + _nodes_side_B)
	r, c = triu(A).nonzero()
	edges = np.array([[rc[0], rc[1]] for rc in zip(r, c)]).astype(int)
	Np = len(_nodes_side_A)
	Nr = len(_nodes_side_B)
	
	# Pass the edge list to a c++ function (src/_km_ompnet.cpp)	
	results = _cp._detect(edges = edges,\
			_nodes_side_A = np.array(list(range(Np))).astype(int),\
			_nodes_side_B = np.array(list(range(Np, Nr + Np))).astype(int),\
			node_capacity = np.array(node_capacity).astype(float),\
			resol = float(resol),\
			num_samples = int(num_samples),\
			num_runs = 10,\
			consensus_threshold = float(consensus_threshold),\
			significance_level = float(significance_level),\
			num_rand_nets = int(num_rand_nets))


	""" RETRIEVE THE RESULTS """
	c = results[0].astype(int)	
	x = results[1]
	
	# Exclude homeless nodes that do not belong to any consensus CP pairs  	
	b = c>=0
	c = dict(zip(compress(_nodes_side_A, b), c[b]))	
	x = dict(zip(compress(_nodes_side_A, b), x[b]))
	
	return c, x	

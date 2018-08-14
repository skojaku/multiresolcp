import _km_omp as _cp
import networkx as nx
from itertools import compress
import numpy as np
import scipy
from scipy.sparse import triu 

def detect(G, ports, resol = 1, phi = {}, num_samples = 100, num_runs=10, consensus_threshold=0.9, significance_level = 0.05, num_rand_nets = 500):
	"""
	Parameters
	----------
	
	G: NetworkX graph
	    Bipartite network composed of N nodes of one type and M nodes of another type.
	    See details in `NetworkX documentation <https://networkx.github.io/documentation/stable/>`_.
	
	ports: list of length N
	    Nodes to project (e.g., specify port nodes to create a network of ports)
	
	resol : float (Optional; Default = 1; 0<=resol)
	    Resolution parameter 
	
	phi : dict of length M (Optional; Default phi[route] = 1 for all routes)
	    - key : route name
	    - value : container capacity 
	
	num_samples: int (Optional; Default = 100; 0 < num_samples)
	    Number of sample CP structures used to obtain consensus CP structure
	
	num_runs: int (Optional; Default = 10; 0 < num_runs)
	    Number of runs of the algorithm to find one sample CP structure
	
	consensus_threshold: float (Optional; Default = 0.9; 0 <= consensus_threshold <=1)
	    Consensus threshold
	
	significance_level: float (Optional; Default = 0.05; 0 < significance_level <=1)
	    Statistical significance level before the Šidák correction
	
	num_rand_nets: int (Optional; Default = 500; 0 < num_rand_nets)
	    Number of randomised networks used to infer the statistical significance
	
	Returns
	-------
	
	c: dict of length N
	    - key: port name
	    - value: index of the consensus CP pair to which the port belongs  
	
	x: dict of length N
	    - key: port name
	    - value: coreness of the port
	"""

	# Make a list of route nodes.  	
	routes = [node for node in G.nodes() if node not in ports]
	
	Np = len(ports) # Number of ports
	Nr = len(routes) # Number of routes
	
	# Convert phi (dict) to numpy array. 
	if len(phi) == 0:
		phi = np.array(np.ones(len(routes))).astype(float)
	else:
		phi = np.array([phi[r] for r in routes]).astype(float)

	# Make the list of edges in the given network
	A = nx.adjacency_matrix(G, ports+routes)
	r, c = triu(A).nonzero()
	edges = np.array([[rc[0], rc[1]] for rc in zip(r, c)]).astype(int)
	
	# Pass parameters to a c++ function (src/_km_ompnet.cpp)	
	results = _cp._detect(edges = edges,\
			ports = np.array(list(range(Np))).astype(int),\
			routes = np.array(list(range(Np, Nr + Np))).astype(int),\
			phi = np.array(phi).astype(float),\
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
	c = dict(zip(compress(ports,b), c[b]))	
	x = dict(zip(compress(ports,b), x[b]))
	
	return c,x	

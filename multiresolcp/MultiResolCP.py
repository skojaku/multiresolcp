import _km_ompnet as _cp
import networkx as nx
from itertools import compress
import numpy as np
import scipy
from scipy.sparse import triu 

def detect(G, ports, resol = 1, phi = {}, num_samples = 100, num_runs=10, consensus_threshold=0.9, significance_level = 0.05, num_rand_nets = 500):
	""" Detect core-periphery pairs at each resolution parameter 

	Parameters
	----------
	G : NetworkX graph
	    Bipartite network composed of N nodes of one part and M nodes of the other part.
	
	ports : list of length N
		Nodes to project onto (e.g., specify port nodes to create a network of ports)
	
	resol : float (Optional. Default = 1)
		Resolution parameter  

	phi : dict of length M (Optional. Default phi[route] = 1 for all routes)
		- key : Route name
		- value : Container capacity 
	
	num_samples: float (Optional. Default = 0.9)
		Number of results used to obtain consensus CP structure
	
	num_runs: float (Optional. Default = 0.9)
		Number of runs of the algorithm to obtain one result
	
	consensus_threshold: float (Optional. Default = 0.9)
		Conensus threshold. Range [0,1].
	
	significance_level: float (Optional. Default = 0.05)
		Statistical significance level before the the Šidák correction.
	
	num_rand_nets: float (Optional. Default = 500)
		Number of randomised networks to infer the statistical significance
	
	Returns
	-------
	c : dict of length N
		- key : port name
		- value : index of the CP pair to which the port belongs.  

	x : dict of length N
		- key : port name
		- value : coreness of port.
	"""
	
	routes = [node for node in G.nodes() if node not in ports]

	if len(phi) == 0:
		phi = np.array(np.ones(len(routes))).astype(float)
	else:
		phi = [phi[r] for r in routes]

	A = nx.adjacency_matrix(G, ports+routes)
	r, c = triu(A).nonzero()
	edges = np.array([ [rc[0], rc[1]] for rc in zip(r, c)]).astype(int)
	
	Np = len(ports)
	Nr = len(routes)

	results = _cp._detect(edges = edges, ports = np.array(list(range(Np))).astype(int), routes = np.array(list(range(Np, Nr + Np))).astype(int), phi = np.array(phi).astype(float), resol = float(resol), num_samples = int(num_samples), num_runs = int(num_runs), consensus_threshold = float(consensus_threshold), significance_level = float(significance_level),num_rand_nets = int(num_rand_nets))

	c = results[0].astype(int)	
	x = results[1]	
	b = c>=0
	c = dict(zip(compress(ports,b), c[b]))	
	x = dict(zip(compress(ports,b), x[b]))
	
	return c,x	

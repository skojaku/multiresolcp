import _km_ompnet as _cp
import networkx as nx
import numpy as np
from scipy.sparse import diags
from networkx.algorithms import bipartite
import scipy

class KM_multiresol():

	""" Initialiser """	
	def __init__(self):
		self.num_runs = 10
		self.num_samples = 100
		self.th = 0.90
		self.significance_level = 0.05
		self.x_ = {} 
		self.c_ = {}
		self.Q_ = {}
		self.qs_ = {}
	
	
	def detect(self, G, routes, resols = [1], phi = {}):
		""" Detect core-periphery pairs at each resolution parameter 
	
		Parameters
		----------
		G : NetworkX graph
		    Bipartite network composed of N nodes of one part and M nodes of the other part.
		
		ports : list of length N
			Nodes to project onto (e.g., specify port nodes to create a network of ports)
		
		ports : list of length M
			Nodes to be aggregated (e.g., specify route nodes to create a network of ports)

		resols : list 
			Resolution parameters 
	
		phi : dict of length M (optional)
			phi takes a node name as key and container capacity as its value, e.g., phi[r] is the container capacity of route r.
			If not given, then phi[r] = 1 for all r.
		"""
		
		if len(phi) == 0:
			phi = np.ones(len(routes))
		else:
			phi = [phi[r] for r in routes]

		ports = [node for node in G.nodes() if node not in routes]	
		A = nx.adjacency_matrix(G, ports+routes)
		r, c = A.nonzero()
		edges = [ [rc[0], rc[1]] for rc in zip(r, c)]
	
		Np = len(ports)
		Nr = len(routes)
		
		# Detect CP pairs at each resolution parameter
		clist = []
		xlist = []
		for resol in resols:
			print(resol)
			results = _cp._detect(edges, np.array(range(Np)).astype(int), np.array(range(Np, Nr + Np)).astype(int), np.array(phi).astype(float), resol )
			
			print(results)	
				
		self.c_ = dict(zip(resols, clist))
		self.x_ = dict(zip(resols, xlist))
		
		
		
	def _matching(self, clist):
		return 	
			
		
		
	
	""" Setters """	

	def set_resolutions(self, resols):
		self.resols = resols
		
	
	""" Getters """

	def get_pair_id(self):
		"""Get core-periphery pair ID of each node.
	
		
		Returns
		-------
		c : dict
			Key: Node name
			Value: IDs of core-periphery pair to which it belongs. 
			
		"""
		return self.c_

	def get_coreness(self):
		"""Get core-periphery pair ID of each node.
	
		
		Returns
		-------
		x : dict
			Key: Node name
			Value: IDs of core-periphery pair to which it belongs. 
			
		"""
		return self.x_
	def score(*args):
		"""Get score of core-periphery pairs.
		
		Parameters
		----------
		G : Graph object. 
		c : Dict object,  the keys and values of which are the name of node and its ID of belonging core-periphery pair.  
		
		
		Returns
		-------
		q : List. q[i] is the quality of core-periphery pair i.  
			
		"""
		self = args[0]

		if len(args) ==1:
			return self.qs_
		else:
			G = args[1]
			c = args[2]
			x = args[3]
			return self._score(G, c, x)	
		

	def _score(self, G, c, x):

		node_pairs, w, node2id, id2node = _to_edge_list(G)
	
		N = len(id2node)
		_c = np.array([ c[id2node[i]]  for i in range(N) ])
		_x = np.array([ x[id2node[i]]  for i in range(N) ])
			
		result = _cp.calc_Q_modmat(edges=node_pairs, ws=w, c = _c, x=_x)

		return result[1].tolist()

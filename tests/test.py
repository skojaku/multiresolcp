import networkx as nx
import numpy as np
import pandas as pd
import multiresolcp as mcp 


def test():
	# Read edge list (space-separated file)
	df = pd.read_csv('example/data/edge-list.dat', sep=' ')
	
	# Read the capacity of each route 
	df2 = pd.read_csv('example/data/capacity.dat', sep=' ')
	
	# Construct NetworkX graph object
	G = nx.from_pandas_edgelist(df)
	
	# Make a dict object of capacities 
	capacity = dict(zip(df2.name.values, df2.capacity.values))
	
	# Make a list of port nodes 
	ports = df['source'].unique().tolist()
	
	# Detect core-periphery structure of the network of ports.
	c, x = mcp.detect(G, ports, resol = 1, phi = capacity, consensus_threshold = 0.9)
	
	# Show the detected consensus CP pairs 
	for k in list(c.keys()):
		print('%s: %d %f' % (k, c[k], x[k]))

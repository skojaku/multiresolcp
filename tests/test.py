import networkx as nx
import numpy as np
import pandas as pd
import multiresolcp as mcp 


def test():
	
	df = pd.read_csv('example/data/edge-list.dat', sep=' ') # Load a list of edges (space-separated file)
	
	G = nx.from_pandas_edgelist(df) # NetworkX graph object
	
	part1 = df['source'].unique().tolist() # List of nodes in part 1
	part2 = df['target'].unique().tolist() # List of nodes in part 2
	
	for resolution in [0.01, 0.1, 0.5, 1, 1.5, 2]:	
	
		c, x = mcp.detect(G, part1, part2, part_to_project = 'part1', resol = resolution) # Detect core-periphery structure at 'resolution'
	
		# Show results	
		print('')
		print('resolution = %f' % resolution)
		print('c:', c)
		print('x:', x)

import networkx as nx
import numpy as np
import pandas as pd
import multiresolcp as mcp 


def test():
	filename = 'samples/unicodelang/unicode.dat' 
	
	df = pd.read_csv(filename, sep='\t', index_col = False, header=None, names = ['source', 'target'], keep_default_na=False, na_values=[''])
	
	ports = df['source'].unique().tolist()
	routes = df['target'].unique().tolist()
	
	G = nx.from_pandas_edgelist(df)
	c, x = mcp.detect(G, routes, resol = 2, consensus_threshold = 0.9)
	
	for k in list(c.keys()):
		print('%s: %d %f' % (k, c[k], x[k]))

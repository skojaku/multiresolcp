import networkx as nx
import numpy as np
import pandas as pd
import multiresolcp as mcp 

filename = 'samples/unicodelang/unicode.dat' 
#filename = 'samples/brunson_revolution/out.brunson_revolution_revolution' 
#df["target"] = df['target']+df['source'].max()

df = pd.read_csv(filename, sep='\t', index_col = False, header=None, names = ['source', 'target'], keep_default_na=False, na_values=[''])

ports = df['source'].unique().tolist()
routes = df['target'].unique().tolist()

G = nx.from_pandas_edgelist(df)
print(len(G.edges()))
c, x = mcp.detect(G, routes, resol = 2, consensus_threshold = 0.9)
print(c)

for k in list(c.keys()):
	print('%s: %d %f' % (k, c[k], x[k]))

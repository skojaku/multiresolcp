import networkx as nx
import numpy as np
import pandas as pd
import multiresolcp as mcp 

df = pd.read_csv('samples/sample.dat', sep='\t', index_col = False)
ports = df["source"].unique().tolist()
routes = df["target"].unique().tolist()
G = nx.from_pandas_edgelist(df)

c, x = mcp.detect(G, routes, resol = 2)

for k in list(c.keys()):
	print('%s: %d %f' % (k, c[k], x[k]))

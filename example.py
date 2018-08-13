import networkx as nx
import numpy as np
import pandas as pd
import multiresolcp as mcp 

df = pd.read_csv('test/data/port_vs_route_w1.dat', sep='\t', index_col = False)
df["source"]-=1
df["target"]-=1
df["target"]+=df["source"].max() + 1
ports = df["source"].unique().tolist()
routes = df["target"].unique().tolist()
G = nx.from_pandas_edgelist(df)

df = pd.read_csv('test/data/route_stat.dat', sep='\t', index_col = False)
phi = df["capacity"].tolist()
phi = dict(zip(routes, phi))
 
c, x = mcp.detect(G, ports, resol = 1)


#print(c, x)

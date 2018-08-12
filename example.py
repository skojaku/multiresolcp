import networkx as nx
import numpy as np
import pandas as pd
import km_ompnet

df = pd.read_csv('samples/out.brunson_club-membership_club-membership', sep='\t', index_col = False)
df["source"]-=1
df["target"]-=1
df["target"]+=df["source"].max() + 1

df2 = df.copy()
df2["source"]+=40
df2["target"]+=40
df = pd.concat([df, df2], ignore_index=True)

G = nx.from_pandas_edgelist(df)

km = km_ompnet.KM_multiresol()

ports = list(range(25))
routes = list(range(25, 40))

ports = ports + [p + 40 for p in ports] 
routes = routes + [r + 40 for r in routes] 

phi = (np.ones(len(routes))).tolist()
phi = dict(zip(routes, phi))

#np.set_printoptions(threshold=np.nan)
km.detect(G, routes, resols = list(range(1,10)))



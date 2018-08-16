import networkx as nx
import numpy as np
import pandas as pd
import multiresolcp as mcp 

# Read edge list (space-separated file)
df = pd.read_csv('data/edge-list.dat', sep=' ')

# Construct NetworkX graph object
G = nx.from_pandas_edgelist(df)

# Make a list of nodes in part 1 and that of nodes in part 2
part1 = df['source'].unique().tolist()
part2 = df['target'].unique().tolist()

# Make a dict object of node capacities 
df2 = pd.read_csv('data/capacity.dat', sep=' ')
node_capacity = dict(zip(df2.name.values, df2.capacity.values))

# Detect core-periphery structure in the network of nodes in part 1 
c, x = mcp.detect(G, part1, part2, part_to_project = 'part1', node_capacity = node_capacity)

# Show the detected consensus CP pairs 
print('Core-periphery structure in the network of nodes in part 1')
for k in  sorted(c, key=c.get):
	print('%s: %d %f' % (k, c[k], x[k]))
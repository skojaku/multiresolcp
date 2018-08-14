import csv
import numpy as np
import pandas as pd
import networkx as nx
import _cpalgorithm as _cp
import cpalgorithm as cp

G=nx.karate_club_graph()
#G=nx.florentine_families_graph()
#df = pd.read_csv("karate.dat", sep='\t');
#G = nx.from_pandas_edgelist(df, "source", 'target', 'weight')

be = cp.BE()

Q = []
be.detect(G)
c = be.get_pair_id()
x = be.is_core()

print(sum(be.score()))

#significance, p_values, q_tilde, s_tilde = cp.qstest(c, x, G, be, num_of_thread = 4, null_model = cp.erdos_renyi)
print(c,x)
#print(significance, p_values)

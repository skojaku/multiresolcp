
About this library
==================

A Python wrapper for the c++ code of the KM algorithm for networks induced by a one-mode projection of bipartite networks.

Please cite:

Sadamori Kojaku, Mengqiao Xu, Haoxiang Xia and Naoki Masuda.
Multiscale core-periphery structure in a global liner shipping network.
Preprint arXiv: ????.

Contents
========
Python code:
  - multiresolcp/__init__.py
  - multiresolcp/MultiResolCP.py

C++ code: 
  - include/km_omp.h
  - include/km_multiresol.h
  - include/graph.h

Python - C++ interface:
  - src/_km_omp.h

Example data and code:
  - example.py
  - data/edge-list.dat 
  - data/capacity.dat 

Others (for PyPi registration and Travis-CI):
  - MANIFEST.in
  - requirements.txt
  - setup.py
  - .travis.yml
  - tests

Installation
============

To install, type
      
.. code-block:: bash

  pip3 install multiresolcp 

If you don't have root privilege, use -user flag, i.e.,  
      
.. code-block:: bash

  pip3 install --user multiresolcp 


Usage
=====

.. code-block:: bash
  
  import multiresolcp
  c, x = multiresolcp.detect(G, ports, resol, phi, num_samples, num_runs, consensus_threshold, significance_level, num_rand_nets)

Parameters
----------

G: NetworkX graph
    Bipartite network composed of N nodes of one type and M nodes of another type.
    See details in `NetworkX documentation <https://networkx.github.io/documentation/stable/>`_.

ports: list of length N
    Nodes to project (e.g., specify port nodes to create a network of ports)

resol : float (Optional; Default = 1; 0<=resol)
    Resolution parameter 

phi : dict of length M (Optional; Default phi[route] = 1 for all routes)
    - key : route name
    - value : container capacity 

num_samples: int (Optional; Default = 100; 0 < num_samples)
    Number of sample CP structures used to obtain consensus CP structure

num_runs: int (Optional; Default = 10; 0 < num_runs)
    Number of runs of the algorithm to find one sample CP structure

consensus_threshold: float (Optional; Default = 0.9; 0 <= consensus_threshold <=1)
    Consensus threshold

significance_level: float (Optional; Default = 0.05; 0 < significance_level <=1)
    Statistical significance level before the Šidák correction

num_rand_nets: int (Optional; Default = 500; 0 < num_rand_nets)
    Number of randomised networks used to infer the statistical significance

Returns
-------

c: dict of length N
    - key: port name
    - value: index of the consensus CP pair to which the port belongs  

x: dict of length N
    - key: port name
    - value: coreness of the port

Note that c and x only contain the nodes in the consensus CP pairs.
If c and x do not contain some nodes, it means that these missing nodes do not belong to any consensus CP pair. 
If you obtain too few nodes in c and x, try decreasing the consensus threshold (i.e., consensus_threshold).
    

Example
=======

.. code-block:: python
  
  import networkx as nx
  import numpy as np
  import pandas as pd
  import multiresolcp as mcp 
  
  # Read edge list (space-separated file)
  df = pd.read_csv('data/edge-list.dat', sep=' ')
  
  # Read the capacity of each route 
  df2 = pd.read_csv('data/capacity.dat', sep=' ')
  
  # Construct NetworkX graph object
  G = nx.from_pandas_edgelist(df)
  
  # Make a dict object of capacities 
  capacity = dict(zip(df2.name.values, df2.capacity.values))
  
  # Make a list of port nodes 
  ports = df['source'].unique().tolist()
  
  # Detect core-periphery structure of the network of ports.
  c, x = mcp.detect(G, ports, resol = 1, phi = capacity, consensus_threshold = 0.9, significance_level = 1.0)
  
  # Show the detected consensus CP pairs 
  for k in list(c.keys()):
  	print('%s: %d %f' % (k, c[k], x[k]))

Requirements
============
- Python 3.4 or later
- Numpy 1.14 or later
- SciPy 1.1 or later
- NetworkX 2.0 or later
- pybind11 2.2 or later 

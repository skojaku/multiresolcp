
About this library
==================

Python code of the KM algorithm for networks induced by a one-mode projection of bipartite networks.

Please cite:

Multiscale core-periphery structure in a global liner shipping network
Sadamori Kojaku, Mengqiao Xu, Haoxiang Xia and Naoki Masuda
Preprint arXiv: ????

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

resol : float (Optional. Default = 1)
    Resolution parameter 

phi : dict of length M (Optional. Default phi[route] = 1 for all routes)
    - key : route name
    - value : container capacity 

num_samples: int (Optional. Default = 100. 0 < num_samples )
    Number of sample CP structures used to obtain consensus CP structure

num_runs: int (Optional. Default = 10. 0 < num_runs )
    Number of runs of the algorithm to find one sample CP structure

consensus_threshold: float (Optional. Default = 0.9. 0 <= consensus_threshold <=1 )
    Consensus threshold

significance_level: float (Optional. Default = 0.05. 0 < significance_level <=1 )
    Statistical significance level before the Šidák correction

num_rand_nets: int (Optional. Default = 500. 0 < num_rand_nets )
    Number of randomised networks used to infer the statistical significance

Returns
-------

c: dict of length N
    - key: port name
    - value: index of the consensus CP pair to which the port belongs.  

x: dict of length N
    - key: port name
    - value: coreness of the port.

Example
=======

.. code-block:: python
  
  import multiresolcp
  import networkx as nx
  
  # Load NetworkX graph object from file "test.edgelist"
  G = nx.read_edgelist("test.edgelist")

  # Nodes to project (e.g., port nodes if one intends to create a network of ports) 
  ports = [0, 1, 2, 3, 4] 

  # Container capacity (e.g., key = route name, value = container capacity) 
  phi = {5 : 0.5, 6 : 0.1, 7 : 0.8, 8: 0.2, 9 : 1} 
  
  # Detect core-periphery structure  
  c, x = multiresolcp.detect(G, ports, phi = phi)

  # Show results 
  print('Pair\tCoreness')
  for k in c.keys():
    print('%d\t%f' % (c[k], x[k]))

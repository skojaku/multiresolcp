
About this library
==================

Python code for the KM algorithm for networks constructed from bipartite networks.

Please cite

Multiscale core-periphery structure in a global liner shipping network
Sadamori Kojaku, Mengqiao Xu, Haoxiang Xia and Naoki Masuda
Preprint arXiv: ????

Installation
============

To install, type
      
.. code-block:: bash

  pip3 install multiresolcp 

If you don't have a root privilege, use -user flag, i.e.,  
      
.. code-block:: bash

  pip3 install --user multiresolcp 

If you don't have a root privilege, use -user flag, i.e.,  
      
.. code-block:: bash

  pip3 install --user multiresolcp 

Usage
=====

.. code-block:: bash
  
  import multiresolcp
  c, x = multiresolcp.detect(G, ports, resol, phi, num_results, num_runs, consensus_threshold, significance_level, num_rand_nets)

Parameters
----------

G : NetworkX graph
    Bipartite network composed of N nodes of one part and M nodes of the other part.
    See details in `NetworkX documentation <https://networkx.github.io/documentation/stable/>`_.

ports : list of length N
	Nodes to project onto (e.g., specify port nodes to create a network of ports)

resol : float (Optional. Default = 1)
	Resolution parameter  

phi : dict of length M (Optional. Default phi[route] = 1 for all routes)
	- key : Route name
	- value : Container capacity 

num_results: float (Optional. Default = 100)
	Number of results used to find consensus CP structure

num_runs: float (Optional. Default = 10)
	Number of runs of the algorithm for one result

consensus_threshold: float (Optional. Default = 0.9)
	Conensus threshold. Range [0,1].

significance_level: float (Optional. Default = 0.05)
	Statistical significance level before the Šidák correction. Range [0,1]

num_rand_nets: float (Optional. Default = 500)
	Number of randomised networks to infer the statistical significance

Returns
-------

c : dict of length N
	- key : port name
	- value : index of the CP pair to which the port belongs.  

x : dict of length N
	- key : port name
	- value : coreness of port.

Example
=======

.. code-block:: bash
  
  import multiresolcp 
  Write example code here
  c, x = multiresolcp.detect(G, ports, resol, phi, num_results, num_runs, consensus_threshold, significance_level, num_rand_nets)

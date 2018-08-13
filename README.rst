Detect core-periphery pairs at each resolution parameter 

Parameters
----------
G : NetworkX graph
    Bipartite network composed of N nodes of one part and M nodes of the other part.

ports : list of length N
	Nodes to project onto (e.g., specify port nodes to create a network of ports)

resol : float (Optional. Default = 1)
	Resolution parameter  

phi : dict of length M (Optional. Default phi[route] = 1 for all routes)
	- key : Route name
	- value : Container capacity 

num_results: float (Optional. Default = 0.9)
	Number of results used to obtain consensus CP structure

num_runs: float (Optional. Default = 0.9)
	Number of runs of the algorithm to obtain one result

consensus_threshold: float (Optional. Default = 0.9)
	Conensus threshold. Range [0,1].

significance_level: float (Optional. Default = 0.05)
	Statistical significance level before the the Šidák correction.

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

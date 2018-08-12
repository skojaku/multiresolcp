import numpy as np
import networkx as nx
import multiprocessing as mp
from scipy.stats import norm

def sz_n(network, c, x):
    return np.bincount(list(c.values())).tolist()

def sz_degree(network, c, x):
    K = max(c.values())+1
    w = [0 for i in range(K)]
    for key, val in c.items():
       w[val]+=network.degree(key)

    return w 

def config_model(G):
    deg = [d[1] for d in G.degree()]
    return nx.configuration_model(deg)	

def erdos_renyi(G):
    n = G.number_of_nodes()
    p = nx.density(G)
    return nx.fast_gnp_random_graph(n, p)

def qstest(pair_id, coreness, G, cpa, significance_level=0.05, null_model = config_model, sfunc = sz_n, num_of_thread = 4, num_of_rand_net = 500 ):
    """(q,s)-test for core-periphery structure.
    
    This function computes the significance of individual core-periphery pairs using either the Erdos-Renyi or the configuration model as the null model. 
    
    Parameters
    ----------
    pair_id : dict
	keys and values of which are node names and IDs of core-periphery pairs, respectively.
    
    coreness : dict
	keys and values of which are node names and coreness, respectively. 
    
    G : NetworkX graph object 
	
    cpa : CPAlgorithm class object
	Core-periphery detection algorithm 
    
    significance_level : float
	Significance level (optional, default 0.5) 
    
    null_model : function
	Null model for generating randomised networks.
       	Provide either config_model or erdos_renyi (optional, default config_model). 
       	One can use another null models. 
       	Specifically, one needs to define a function taking NetworkX graph object as input and randomised network as its output. 
       	Then, one gives the defined function, say myfunc,  to qstest by null_model=myfunc.
    
    sfunc : function
	Size function (optional, default sz_n)
       In the (q,s)--test, one is required to provide a function for measuring the size of an individual core-periphery pair. By default, this function is the number of nodes in the core-periphery pair (i.e., sz_n). One can set sz_degree, which measures the size as the sum of the degree of nodes belonging to the core-periphery pair.  
    
    num_of_thread : function
	Number of thread (optional, default 4)
     
    	The (q,s)--test uses multiple threads to compute the significance. 
    
    num_of_rand_net : int
	Number of randomised networks (optional, default 500)
    
    Returns
    -------
    sig_pair_id : dict
	keys and values of which are node names and IDs of core-periphery pairs, respectively. If nodes belong to insignificant core-periphery pair, then the values are None. 

    sig_coreness : dict 
	significance[i] = True or significance[i] = False indicates core-periphery pair i is significant or insignificant, respectively. If nodes belong to insignificant core-periphery pair, then the values are None.

    significance : list 
	significance[i] = True or significance[i] = False indicates core-periphery pair i is significant or insignificant, respectively.
    
    p_values : list
	p_values[i] is the p-value of core-periphery pair i.
    
    Examples
    --------
    Detect core-periphery pairs in the karate club network.
    
    >>> import cpalgorithm as cpa	
    >>> km = cpa.KM_config()
    >>> km.detect(G) 
    >>> pair_id = km.get_pair_id() 
    >>> coreness = km.get_coreness()
    
    Examine the significance of each core-periphery pair using the configuration model:	
    
    >>> sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, km)
    
    or
    
    >>> sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, km, null_model=config_model)
    
    Examine the significance of each core-periphery pair using the Erdos-Renyi random graph:
    
    >>>  sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, km, null_model=erdos_renyi)
    	
    .. rubric:: Reference
    
    Sadamori Kojaku and Naoki Masuda.
    A generalised significance test for individual communities in networks.
    Scientific Reports, 8:7351 (2018)
    """
   
    q = np.array(cpa.score(G, pair_id, coreness), dtype = np.float)    
    s = np.array(sfunc(G, pair_id, coreness) , dtype = np.float)
    C = len(q)
    alpha_corrected = 1.0 - (1.0 - significance_level) ** (1.0 / float(C))
        
    q_tilde = []
    s_tilde = []
    if num_of_thread == 1:
        q_tilde, s_tilde = draw_qs_samples(G, sfunc, cpa, null_model, num_of_rand_net)
    else:
        private_args = [(G, sfunc, cpa, null_model, int(num_of_rand_net / num_of_thread) + 1) for i in range(num_of_thread)]
        pool = mp.Pool(num_of_thread)
        qs_tilde = pool.map(wrapper_draw_qs_samples, private_args)
        for i in range(num_of_thread):
            q_tilde += qs_tilde[i][0] 
            s_tilde += qs_tilde[i][1]
    
    q_tilde = np.array(q_tilde, dtype = np.float)    
    s_tilde = np.array(s_tilde, dtype = np.float)    
    q_ave = np.mean(q_tilde)
    s_ave = np.mean(s_tilde)
    q_std = np.std(q_tilde, ddof = 1)
    s_std = np.std(s_tilde, ddof = 1)
   
    if (s_std <= 1e-30) or (q_std <= 1e-30):
        gamma = 0.0
        s_std = 1e-20
    else:
        gamma = np.corrcoef(q_tilde, s_tilde)[0, 1]
     
    h = float(len(q_tilde)) ** (- 1.0 / 6.0)
    p_values = [1.0] * C
    significant = [False] * C

    cidx = 0
    cid2newcid = - np.ones(C)
    for cid in range(C):
        if (s_std <= 1e-30) or (q_std <= 1e-30):
            continue    
        w = np.exp(- ( (s[cid] - s_tilde) / (np.sqrt(2.0) * h * s_std) ) ** 2)
        cd = norm.cdf( ( (q[cid] - q_tilde) / (h * q_std) - gamma * (s[cid] - s_tilde) / (h * s_std) ) / np.sqrt(1.0 - gamma * gamma) )    
        denom = sum(w)    
        if denom <= 1e-30:
            continue    
        p_values[cid] = 1.0 - (sum( w * cd ) / denom)
        significant[cid] = p_values[cid] <= alpha_corrected
        
        if significant[cid]:
            cid2newcid[cid] = cidx
            cidx+=1 
        
    sig_pair_id = pair_id
    sig_coreness = coreness
   
	 
    for k, v in sig_pair_id.items():
        if significant[v]:
            sig_pair_id[k]=cid2newcid[ pair_id[k] ]
        else:
            sig_pair_id[k]=None
            sig_coreness[k]=None
        
    return sig_pair_id, sig_coreness, significant, p_values


# Private function for qstest        
def draw_qs_samples(G, sfunc, cpa, null_model, num_of_rand_net):
    #deg = [x[1] for x in G.degree()]
    q_rand = []
    s_rand = []

    for i in range(num_of_rand_net):
        Gr = null_model(G)
        cpa.detect(Gr) 
        q_rand = q_rand + cpa.score()
        s_rand = s_rand + sfunc(Gr, cpa.get_pair_id(), cpa.get_coreness()) 
    return q_rand, s_rand


# Private function for qstest        
def wrapper_draw_qs_samples(args):
    return draw_qs_samples(*args)    

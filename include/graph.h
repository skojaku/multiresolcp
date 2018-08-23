/*

This code defines two class objects, AdjacentNode and Graph.

AdjacentNode is an alias of c++ pair class.  

Graph contains information on the weighted adjacency matrix of an undirected network.
Graph also accommodates functions for computing some basic statistics such as the number of nodes, degree of each node and weighted degree of each node.

*/

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

class AdjacentNode {
 public:
  int node;
  double weight;
  AdjacentNode(int _node, double _weight) {
    node = _node;
    weight = _weight;
  }
};


class Graph {
 public:

  /* Constructor */
  Graph() { _N = 0; }

  /* 
  
  Getters
  ------
  get_num_nodes() : Return the number of nodes in the network. 

  get_total_edge_weight() : Return the sum of the weight of edges in the network.  

  degree(int nid) : Return the degree of node nid.  

  wdegree(int nid) : Return the sum of the weight of edges emanating from node nid 
	
  contains(int nid) : Return true if the network contains the node with ID nid. Otherwise, return false.  

  adjacency_list() : Return the list of neighbours of each node. 
                     Specifically, 'adjacency_list()[i]' is the list of neighbours of node i.
                     The 'adjacency_list()[i]' is a vector object whose jth element is an AdjacentNode object containing 
	             the ID of the neighbouring node (i.e., 'adjacency_list()[i][j].node') and the weight of the edge 
                     between node i and the neighbouring node (i.e., 'adjacency_list()[i][j].weight').

  neighbours(int nid) : Return the list of neighbours of node i. 
                        The jth element is an AdjacentNode object containing the ID of the neighbouring node 
                        (i.e., 'adjacency_list()[i][j].node') and the weight of the edge between node i and the neighbouring 
                        node (i.e., 'adjacency_list()[i][j].weight').

  */
  int get_num_nodes();
  int get_num_edges(); 
  double get_total_edge_weight();
  int degree(int nid);
  double wdegree(int nid);
  bool contains(int nid);
  map<int, vector<AdjacentNode>>& adjacency_list();
  vector<AdjacentNode>& neighbours(int nid);

  /* 
   
  Add an undirected edge between nodes u and v. 
 
 
  Parameters 
  ----------
  u : ID of node  

  v : ID of a neighbouring node of node u (swapping u and v does not change anything)
 
  w : Weight of the edge 

  */
  void addEdge(int u, int v, double w);

  /* 

  Replace multiple edges between the same pair of nodes with a single edge whose weight is 
  is the sum of the weight of the multiple edges.
  
  */
  void aggregate_multi_edges();

  /* 

  This function replaces the IDs of nodes with new consecutive IDs starting from 0. 
  
  */
  map<int, int> renumbering();

  /* 

  Remove all the nodes and edges in the network. 
  
  */
  void clear();


 /* Private variables */
  private:
   map<int, vector<AdjacentNode>> _A;
   int _N;
   vector<AdjacentNode> _empty_vec;
};

int Graph::get_num_nodes() { return _N; };

int Graph::get_num_edges() {
  int M = 0;
  for (auto& node : _A) {
    M += node.second.size();
  };
  return M / 2;
};

double Graph::get_total_edge_weight() {
  double M = 0;
  for (auto& node : _A) {
    M += wdegree(node.first);
  };
  return M / 2;
};

int Graph::degree(int nid) {
  if (contains(nid)) {
    return (int)_A[nid].size();
  } else {
    return 0;
  }
}

double Graph::wdegree(int nid) {
  if (contains(nid)) {
    double wdeg = 0;
    for (auto& adj : _A[nid]) {
      wdeg += adj.weight;
    };
    return wdeg;
  } else {
    return 0;
  }
}

bool Graph::contains(int nid) {
  if (_A.count(nid) == 0) {
    return false;
  }
  return true;
};

map<int, vector<AdjacentNode>>& Graph::adjacency_list() { return _A; };

vector<AdjacentNode>& Graph::neighbours(int nid) {
  if (contains(nid)) return _A[nid];
  return _empty_vec;
};

void Graph::addEdge(int u, int v, double w) {
  if (0 == _A.count(u)) {
    _A.insert(make_pair(u, vector<AdjacentNode>()));
  }

  if (0 == _A.count(v)) {
    _A.insert(make_pair(v, vector<AdjacentNode>()));
  }

  AdjacentNode adj1(v, w);
  _A[u].push_back(adj1);

  AdjacentNode adj2(u, w);
  _A[v].push_back(adj2);

  _N = max(_N, max(u + 1, v + 1));
}

void Graph::aggregate_multi_edges() {
  Graph Gnew;

  for (auto& node : _A) {
    map<int, double> myMap;
    for (auto& adj : node.second) {
      if (myMap.count(adj.node) == 0) {
        myMap.insert(make_pair(adj.node, adj.weight));
      } else {
        myMap[adj.node] += adj.weight;
      }
    }
    for (const auto& p : myMap) {
      if (node.first > p.first) continue;
      if (node.first == p.first) {
        Gnew.addEdge(node.first, p.first, p.second / 2);
      } else {
        Gnew.addEdge(node.first, p.first, p.second);
      }
    }
  }

  clear();

  _A = Gnew._A;
  _N = Gnew._N;
}

map<int, int> Graph::renumbering() {
  Graph Gnew;
  map<int, int> myMap;
  map<int, int> inv_myMap;
  int idx = 0;
  for (auto& node : _A) {
    myMap.insert(make_pair(node.first, idx));
    inv_myMap.insert(make_pair(idx, node.first));
    idx++;
  }

  for (auto& node : _A) {
    for (auto& adj : node.second) {
      if (node.first > adj.node) continue;
      if (node.first == adj.node) {
        Gnew.addEdge(myMap[node.first], myMap[adj.node], adj.weight / 2);
      } else {
        Gnew.addEdge(myMap[node.first], myMap[adj.node], adj.weight);
      }
    }
  }

  clear();

  _A = Gnew._A;
  _N = Gnew._N;
  return inv_myMap;
}

void Graph::clear() {
  for (auto&& p : _A) {
    p.second.clear();
  }
  _A.clear();
}

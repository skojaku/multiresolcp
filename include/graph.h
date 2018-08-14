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
  map<int, vector<AdjacentNode>> _A;
  int _N;
  vector<AdjacentNode> _empty_vec;

  Graph() { _N = 0; }

  // Getter
  int get_num_nodes() { return _N; };

  int get_num_edges() {
    int M = 0;
    for (auto& node : _A) {
      M += node.second.size();
    };
    return M / 2;
  };

  double get_total_edge_weight() {
    double M = 0;
    for (auto& node : _A) {
      M += wdegree(node.first);
    };
    return M / 2;
  };

  int degree(int nid) {
    if (isnode(nid)) {
      return (int)_A[nid].size();
    } else {
      return 0;
    }
  }

  double wdegree(int nid) {
    if (isnode(nid)) {
      double wdeg = 0;
      for (auto& adj : _A[nid]) {
        wdeg += adj.weight;
      };
      return wdeg;
    } else {
      return 0;
    }
  }

  bool isnode(int nid) {
    if (_A.count(nid) == 0) {
      return false;
    }
    return true;
  };

  map<int, vector<AdjacentNode>>& adjacency_list() { return _A; };
  vector<AdjacentNode>& neighbours(int nid) {
    if (isnode(nid)) return _A[nid];
    return _empty_vec;
  };

  // Setter
  void addEdge(int u, int v, double w);

  // Others
  void aggregate_multi_edges();
  map<int, int> renumbering();
  void clear() {
    for (auto&& p : _A) {
      p.second.clear();
    }
    _A.clear();
  }
};

void Graph::addEdge(int u, int v, double w) {
  if (0 == _A.count(u)) {  // if node u exists
    _A.insert(make_pair(u, vector<AdjacentNode>()));
  }

  if (0 == _A.count(v)) {  // if node u exists
    _A.insert(make_pair(v, vector<AdjacentNode>()));
  }

  AdjacentNode adj1(v, w);
  _A[u].push_back(adj1);

  AdjacentNode adj2(u, w);
  _A[v].push_back(adj2);

  _N = max(_N, max(u + 1, v + 1));
}

// merge multiple-edges with a single weighted edge
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

// Renumbering nodes
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

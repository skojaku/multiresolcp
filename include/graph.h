#include <iostream>
#include <vector>
#include <map>

using namespace std;

class AdjacentNode{
public:
	int node;
	double weight;
    
	AdjacentNode(int _node, double _weight){	
		node = _node;
		weight = _weight;
	}
};
    
class Graph{
public:
    map<int, vector<AdjacentNode>> _A;
    
    // constructer
    Graph(){}
	
    // Getters
    int get_num_nodes() {return _A.size();};

    int get_num_edges() {int M = 0; for(auto& node : _A){M+=node.second.size();}; return M/2;};

    int degree(int nid) {return (int)_A[nid].size();}

    double wdegree(int nid) {int wdeg = 0; for(auto& adj : _A[nid]){wdeg+=adj.weight;}; return wdeg;}
    
    vector<AdjacentNode> &neighbours(int nid){ return _A[nid];};
    map<int, vector<AdjacentNode>> &adjacency_list(){ return _A;};
	
    // Setters
    void addEdge(int u, int v, double w);

    // Others
    void aggregate_multi_edges();

    void clear(){for (auto&& p : _A) {  p.second.clear(); } _A.clear();}

};

    
void Graph::addEdge(int u, int v, double w){
	
	if( 0 < _A.count(u) ){ // if node u exists
		_A.insert(make_pair(u, vector<AdjacentNode>()));
	}

	if( 0 < _A.count(v) ){ // if node u exists
		_A.insert(make_pair(v, vector<AdjacentNode>()));
	}

	AdjacentNode adj1(v, w);
	_A[u].push_back(adj1);
	
	AdjacentNode adj2(u, w);
	_A[v].push_back(adj2);
}

// merge multiple-edges with a single weighted edge 
void Graph::aggregate_multi_edges(){
	
        Graph Gnew;
	
	for(auto& node : _A){
		map<int, double> myMap;
		for(auto& adj : node.second){	
			if ( !myMap.insert( make_pair( adj.node, adj.weight ) ).second ) {
				myMap[adj.node]+=adj.weight;
			}
		}
		for (const auto & p : myMap) {
			Gnew.addEdge(node.first, p.first, p.second);	
		}
	}

	clear();	
	
	_A = Gnew._A;
}


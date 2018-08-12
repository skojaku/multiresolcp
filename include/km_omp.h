/*
*
* Multiresolution version of the KM algorithm algorithm
*
* Multiscale core-periphery structure in a global liner shipping network
* Sadamori Kojaku, Mengqiao Xu, Haoxiang Xia and Naoki Masuda
* Preprint arXiv: ????
*
* Please do not distribute without contacting the authors.
*
* AUTHOR - Sadamori Kojaku
*
* DATE - 08 Aug 2018
*
*/

#include <omp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <vector>

using namespace std;

class KM_omp {
 public:
  // Constructor
  KM_omp(int num_runs, int num_samples);

  void detect(Graph& biG, vector<int>& ports, vector<int>& routes, map<int, double>& phi, double resol);

  vector<int> get_c() const { return _c; };
  vector<double> get_x() const { return _x; };

 private:
  int _num_runs_KM_multiresol;
  int _num_samples;
  int _num_rand_nets;
  double _significance_level;
  vector<double> _nhat;
  vector<double> _qhat;
  vector<int> _c;
  vector<double> _x;
  double _consensus_th;

  void _one_mode_projection(
      Graph& biG, vector<int>& ports, vector<int>& routes, map<int, double>& phi, Graph& uniG);

  void _est_null_model(Graph& biG, vector<int>& ports, vector<int>& routes, map<int, double>& phi, double resol, int num_of_net) ;

  vector<double> _calc_p_values(vector<double>& q, vector<double>& n);
  vector<int> _connected_components(Graph& U, double th, int N);
};

KM_omp::KM_omp(int num_runs, int num_samples) {
  _num_runs_KM_multiresol = num_runs;
  _num_samples = num_samples;
  _num_rand_nets = 500;
  _significance_level = 0.05;
  _consensus_th = 0.9;
}

void KM_omp::detect(
    Graph& biG, vector<int>& ports, vector<int>& routes, map<int, double>& phi, double resol) {

  int Np = ports.size();
  int Nr = routes.size();
  //vector<int> nodelist(Nr + Np,0);
  //for(int i=0; i < Np; i++) nodelist[i]= ports[i]; 
  //for(int i=0; i < Nr; i++) nodelist[i+Np]= routes[i]; 
  
  //biG.reorder(nodelist);

  vector<double> theta(Np, 0.0);
  int idx = 0;
  for (auto& p : ports) {
    theta[idx] = biG.degree(p);
    idx++;
  }

  double phi_d_r = 0;
  double M = biG.get_num_edges();
  for (auto& r : routes) {
    phi_d_r += phi[r] * biG.degree(r);
  }
  resol *= phi_d_r / (M * (M - 1));
  _est_null_model(biG, ports, routes, phi, resol, _num_rand_nets);

  Graph uniG;
  _one_mode_projection(biG, ports, routes, phi, uniG);

  Graph U;
  vector<double> X(Np, 0.0);
  vector<double> sig_count(Np, 0.0);
  KM_multiresol _km(_num_runs_KM_multiresol);
  for (int sid = 0; sid < _num_samples; sid++) {
    _km.detect(uniG, theta, resol);

    vector<int> cs = _km.get_c();
    vector<double> xs = _km.get_x();
    vector<double> qs = _km.get_q();

    int Ks = qs.size();
    vector<double> ns(Ks, 0.0);
    for (int i = 0; i < Np; i++) {
      ns[cs[i]]++;
    }
    vector<double> pvals = _calc_p_values(qs, ns);
    double alpha = 1.0 - pow(1.0 - _significance_level, (1.0 / (double)Ks));
    for (int i = 0; i < Np; i++) {
      if (pvals[cs[i]] > alpha) continue;

      for (int j = i + 1; j < Np; j++) {
        if (pvals[cs[j]] > alpha) continue;

        if (cs[i] == cs[j]) {
          U.addEdge(i, j, 1);
        }
      }
    X[i]+= xs[i];
    sig_count[i]++;
    }
  }
  for (int i = 0; i < Np; i++) {
	if(sig_count[i]==0) continue;
  	X[i]/=sig_count[i];
  }
  U.aggregate_multi_edges();

  _c = _connected_components(U, _num_samples * _consensus_th, Np);
  _x = X;
}

void KM_omp::_one_mode_projection(
    Graph& G, vector<int>& ports, vector<int>& routes, map<int, double>& phi, Graph& uniG) {
  for (auto& r : routes) {
    vector<AdjacentNode> adj = G.neighbours(r);
    int sz = adj.size();
    for (int i = 0; i < sz; i++) {
      for (int j = i + 1; j < sz; j++) {
        uniG.addEdge(adj[i].node, adj[j].node, phi[r] / (double)(G.degree(r) - 1));
      }
    }
  }

  uniG.aggregate_multi_edges();
}

void KM_omp::_est_null_model(Graph& biG, vector<int>& ports, vector<int>& routes, map<int, double>& phi, double resol, int num_of_net) {

  auto _Chung_Lu_Algorithm_bipartite = [](const vector<double>& rdeg, const vector<double>& cdeg,
                                          const vector<int>& rnodes, const vector<int>& cnodes,
                                          vector<int>& rnodes_ids, vector<int>& cnodes_ids, Graph& biG, bool isunweighted,
                                          mt19937_64& mtrnd) {

    uniform_real_distribution<double> udist(0.0, 1.0);
    int Nr = rdeg.size();
    int Nc = cdeg.size();
    double M = accumulate(rdeg.begin(), rdeg.end(), 0.0);
    M /= 2;
    vector<vector<pair<int, double>>> tmp(Nr, vector<pair<int, double>>(0));
    for (int u = 0; u < Nr; u++) {
      for (int v = 0; v < Nc; v++) {
        double p = min(1.0, rdeg[rnodes[u]] * cdeg[cnodes[v]] / (2.0 * M));
        while (v < Nc && p > 0) {
          if (p != 1) {
            double r = udist(mtrnd);
            v = v + floor(log(r) / log(1 - p));
          }
          if (v < Nc) {
            double q = min(rdeg[rnodes[u]] * cdeg[cnodes[v]] / (2.0 * M), 1.0);
            double w = 1;
            bool addEdge = false;
            if (isunweighted) {
              double r = udist(mtrnd);
              addEdge = r < q / p;
            } else {
              poisson_distribution<int> distribution(q / p);
              w = distribution(mtrnd);
              addEdge = w > 0;
            }
            if (addEdge) {
              biG.addEdge(rnodes_ids[rnodes[u]], cnodes_ids[cnodes[v]], w);
            }
            p = q;
            v = v + 1;
          }
        }
      }
    }
  };

  auto _init_random_number_generator = []() {
    mt19937_64 mtrnd;
    random_device r;
    seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    mtrnd.seed(seed);
    return mtrnd;
  };

  /* Main routines */
  bool isunweighted = true;
  int Nport = ports.size();
  int Nroute = routes.size();
  int idx = 0;
  vector<double> rdeg(Nport, 0.0);
  vector<double> cdeg(Nroute, 0.0);
  for (auto& p : ports) {
    rdeg[idx] = biG.degree(p);
    idx++;
  }
  idx = 0;
  for (auto& r : routes) {
    cdeg[idx] = biG.degree(r);
    idx++;
  }

  vector<int> rdeg_rank(Nport);  // deg_rank[k] is the id of the node with the kth largest degree.
  iota(rdeg_rank.begin(), rdeg_rank.end(), 0);
  sort(rdeg_rank.begin(), rdeg_rank.end(), [&](int x, int y) { return rdeg[x] > rdeg[y]; });
  vector<int> cdeg_rank(Nroute);  // deg_rank[k] is the id of the node with the kth largest degree.
  iota(cdeg_rank.begin(), cdeg_rank.end(), 0);
  sort(cdeg_rank.begin(), cdeg_rank.end(), [&](int x, int y) { return cdeg[x] > cdeg[y]; });

  // create random number generator per each thread
  int numthread = 1;
#pragma omp parallel
  { numthread = omp_get_num_threads(); }
  vector<mt19937_64> mtrnd_list(numthread);
  for (int i = 0; i < numthread; i++) {
    mtrnd_list[i] = _init_random_number_generator();
  }
  vector<double> nhat;
  vector<double> qhat;
#ifdef _OPENMP
#pragma omp parallel for shared(mtrnd_list)
#endif
  for (int it = 0; it < _num_rand_nets; it++) {
    int tid = omp_get_thread_num();
    mt19937_64 mtrnd = mtrnd_list[tid];

    Graph biG_rand;
    Graph uniG_rand;
    KM_multiresol _km(_num_runs_KM_multiresol);
    _Chung_Lu_Algorithm_bipartite(rdeg, cdeg, rdeg_rank, cdeg_rank, ports, routes, biG_rand, isunweighted,
                                  mtrnd);
    
    _one_mode_projection(biG_rand, ports, routes, phi, uniG_rand);

    int myidx = 0; 
    vector<double> theta(Nport,0);
    for(auto& p: ports){
	theta[myidx] = biG_rand.degree(p);
	myidx++;
    } 
    _km.detect(uniG_rand, theta, resol);

    vector<double> qr = _km.get_q();

    int K_rand = qr.size();
    vector<int> nsr(K_rand, 0);
    vector<int> cr = _km.get_c();
    for (int i = 0; i < Nport; i++){
      nsr[cr[i]]++;
    }
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      for (int k = 0; k < K_rand; k++) {
        nhat.push_back((double)nsr[k]);
        qhat.push_back(qr[k]);
      }
    }
  }
  _nhat = nhat;
  _qhat = qhat;
}

vector<double> KM_omp::_calc_p_values(vector<double>& q, vector<double>& n) {
  auto normcdf = [](double x) { return 0.5 + 0.5 * erf(x * M_SQRT1_2); };
  int K = q.size();

  int S = _nhat.size();
  double mu_n = (double)accumulate(_nhat.begin(), _nhat.end(), 0.0) / (double)S;
  double mu_q = (double)accumulate(_qhat.begin(), _qhat.end(), 0.0) / (double)S;
  double sig_nn = 0;
  double sig_qq = 0;
  double sig_nq = 0;
  for (int s = 0; s < S; s++) {
    sig_nn += pow((double)_nhat[s] - mu_n, 2) / (double)(S - 1);
    sig_qq += pow(_qhat[s] - mu_q, 2) / (double)(S - 1);
    sig_nq += ((double)_nhat[s] - mu_n) * (_qhat[s] - mu_q) / (double)(S - 1);
  }

  double h = max(pow((double)S, -1.0 / 6.0), 1e-32);
  vector<double> p_values(K, 1.0);
  for (int k = 0; k < K; k++) {
    double numer = 0.0;
    double denom = 0.0;
    for (int s = 0; s < S; s++) {
      double qbar = _qhat[s] + sig_nq / (sig_nn +1e-32) * (double)(n[k] - _nhat[s]);

      double t = sig_nn * (q[k] - qbar) / (sqrt(sig_nn * sig_qq - sig_nq * sig_nq) * h + 1e-32);
      double cum = normcdf(t);

      double w = exp(-(double)pow(n[k] - _nhat[s], 2) / (2.0 * h * h * sig_nn +1e-32)) + 1e-32;
      numer += cum * w;
      denom += w;
    }
    p_values[k] = 1.0 - numer / denom;
  }
  return p_values;
}

vector<int> KM_omp::_connected_components(Graph& U, double th, int N) {
  vector<int> c(N, -1);
  vector<bool> remain(N, false);

  for (auto& nodes : U.adjacency_list()) {
    for (auto& adj : nodes.second) {
      if (adj.weight > th) {
        remain[nodes.first] = true;
        break;
      }
    }
  }
  int cid = 0;
  while (!all_of(remain.begin(), remain.end(), [](bool v) { return !v; })) {
    vector<int> x(N, 0);
    vector<int> xnew(N, 0);
    for (int i = 0; i < N; i++) {
      if (remain[i] == false) continue;
      x[i] = 1;
      break;
    }

    int n = 1;
    int oldn = 0;

    while (oldn != n) {
      oldn = n;

      fill(xnew.begin(), xnew.end(), 0);
      for (auto& nodes : U.adjacency_list()) {
        for (auto& adj : nodes.second) {
          if (adj.weight < th) continue;
          xnew[nodes.first] += x[nodes.first]+x[adj.node];
        }
      }
      x = xnew;
      n = 0;
      for (int i = 0; i < N; i++) {
        if (x[i] > 0) n++;
      }
    }

    for (int i = 0; i < N; i++) {
      if (x[i] > 0) {
        remain[i] = false;
        c[i] = cid;
      }
    }
    cid++;
  }
  for (int i = 0; i < N; i++) {
    if (c[i] < 0) {
      c[i] = cid;
      cid++;
    }
  }
  auto _relabeling = [](vector<int>& c) {
    int N = c.size();
    std::vector<int> labs;
    for (int i = 0; i < N; i++) {
      int cid = -1;
      int labsize = labs.size();
      for (int j = 0; j < labsize; j++) {
        if (labs[j] == c[i]) {
          cid = j;
          break;
        }
      }

      if (cid < 0) {
        labs.push_back(c[i]);
        cid = labs.size() - 1;
      }
      c[i] = cid;
    }
  };
  _relabeling(c);
  return c;
}

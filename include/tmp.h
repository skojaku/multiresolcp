/*
*
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

#include <graph.h>

using namespace std;

class KM_multiresol {
 public:
  // Constructor
  KM_multiresol(int num_runs, vector<double>& phi, double resol);

  // detect and calc_Q take bipartite networks
  void detect(Graph& G, vector<int>& ports, vector<int>& routes);

  void calc_Q(Graph& G,  // bipartite network
              const vector<int>& c,
              const vector<double>& x,
              vector<int>& ports,
              vector<int>& routes,
              const vector<double>& phi,
              double resol,
              vector<double>& q);

  // Getter
  vector<int> get_c() const { return _c; };
  vector<double> get_x() const { return _x; };
  vector<double> get_pvals() const { return _pvals; };
  vector<double> _calc_Q(Graph& G) {
    vector<double> q;
    calc_Q(G, _c, _x, _ports, _routes, _phi, _resol, q);
    return q;
  };

 private:
  /* For core-periphery detection */
  int _num_runs;
  vector<int> _c;
  vector<double> _x;
  double _Q;
  vector<double> _q;
  vector<double> _phi;
  double _resol;

  vector<int> _ports;
  vector<int> _routes;

  /* For statistical test*/
  int _num_rand_nets;
  vector<double> _nhat;
  vector<double> _qhat;
  vector<double> _pvals;

  /* Random number generator */
  mt19937_64 _mtrnd;

  void _label_switching(Graph& G,
                        const int num_of_runs,
                        vector<int>& c,
                        vector<double>& x,
                        double& Q,
                        vector<double>& q,
                        vector<double>& phi,
                        double resol,
                        mt19937_64& mtrnd);

  void _louvain(Graph& G,
                const int num_of_runs,
                vector<int>& c,
                vector<double>& x,
                double& Q,
                vector<double>& q,
                vector<double>& phi,
                double resol,
                mt19937_64& mtrnd);

  void _est_null_model(Graph& G, vector<int>& ports, vector<int>& routes, int num_of_net);

  /* Quality function for projected networks */
  void _calc_Q_projected(Graph& G,
                         const vector<int>& c,
                         const vector<double>& x,
                         map<int, int>& deg_port,
                         double resol,
                         vector<double>& q);

  void _one_mode_projection(
      Graph& G, vector<int>& ports, vector<int>& routes, vector<double>& phi, Graph& G_projected);
  vector<double> _compute_p_values(Graph& G,
                                   vector<int>& c,
                                   vector<double>& x,
                                   vector<int>& ports,
                                   vector<int>& routes,
                                   vector<double>& phi,
                                   double resol);
};

/*-----------------------------
Constructor
-----------------------------*/
KM_multiresol::KM_multiresol(int num_runs, double resol) {
  mt19937_64 mtrnd;
  random_device r;
  seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
  mtrnd.seed(seed);
  _mtrnd = mtrnd;

  _num_runs = num_runs;

  uniform_real_distribution<double> tmp(0.0, 1.0);

  _num_rand_nets = 500;
  _phi = phi;
  _resol = resol;
};

/*-----------------------------
Public functions
-----------------------------*/
void KM_multiresol::detect(Graph& G, vector<int>& ports, vector<int>& routes, vector<double>& phi) {
  Graph G_projected;

  _ports = ports;
  _routes = routes;

  _one_mode_projection(G, ports, routes, _phi, G_projected);

  map<int, int> _deg_port;  
  for(auto& p : ports) _deg_port.insert(make_pair(p, G.degree(p)));

  _louvain(G_projected, _num_runs, _c, _x, _Q, _q, _deg_port, _resol, _mtrnd);

  if (_nhat.size() == 0) {
    _est_null_model(G, _ports, _routes, _num_rand_nets);
  }

  _pvals = _compute_p_values(G, _c, _x, ports, routes, _phi, _resol);
}

void KM_multiresol::_calc_Q_projected(Graph& G,
                                      const vector<int>& c,
                                      const vector<double>& x,
                         		map<int, int>& deg_port,
                                      double resol,
                                      vector<double>& q) {
  int K = *max_element(c.begin(), c.end()) + 1;
  q.assign(K, 0.0);
  vector<double> Dc(K, 0.0);
  vector<double> Dp(K, 0.0);

  double Omega = 0.0;
  for (auto& node : G.adjacency_list()) {
    int i = node.first;
    for (auto& adj : node.second) {
      int j = adj.node;
      double wj = adj.weight;
      q[c[i]] += wj * !!(c[i] == c[j]) * (x[i] + x[j] - x[i] * x[j]);
      Omega += wj;
    }
    Dc[c[i]] += x[i] * deg_port[i];
    Dp[c[i]] += (1 - x[i]) * deg_port[i];
  }
  Omega /= 2;

  for (int k = 0; k < K; k++) {
    q[k] = (q[k] - resol * (Dc[k] * Dc[k] + 2 * Dc[k] * Dp[k])) / (2 * Omega);
  }
}

void KM_multiresol::calc_Q(Graph& G,
                           const vector<int>& c,
                           const vector<double>& x,
                           vector<int>& ports,
                           vector<int>& routes,
                           const vector<double>& phi,
                           double resol,
                           vector<double>& q) {
  Graph G_projected;
  _one_mode_projection(G, ports, routes, _phi, G_projected);
  _calc_Q_projected(G_projected, c, x, phi, resol, q);
}

/*-----------------------------
Private functions
-----------------------------*/
void KM_multiresol::_label_switching(Graph& G,
                                     const int num_of_runs,
                                     vector<int>& c,
                                     vector<double>& x,
                                     double& Q,
                                     vector<double>& q,
                                     vector<double>& phi,
                                     double resol,
                                     mt19937_64& mtrnd) {
  /* Label switching algorithm */
  auto __label_switching = [](Graph& G, vector<int>& c, vector<double>& x, vector<double>& phi, double resol,
                              mt19937_64& mtrnd) {

    auto _propose_new_label = [](Graph& G, const vector<int>& c, const vector<double>& x, const double deg,
                                 const vector<double>& sum_of_deg_core, const vector<double>& sum_of_deg_peri,
                                 const double resol, const int node_id, int& cprime, double& xprime,
                                 double& dQ, mt19937_64& mtrnd) {

      auto _calc_dQ_conf = [](double d_i_c, double d_i_p, double d_i, double D_c, double D_p, double selfloop,
                              double x, const double resol) {
        return 2 * (d_i_c + d_i_p * (x)-resol * d_i * (D_c + D_p * x)) + x * (selfloop - resol * d_i * d_i);
      };

      int N = G.get_num_nodes();

      uniform_real_distribution<double> _udist(0.0, 1.0);

      vector<double> edges_to_core(N, 0.0);
      vector<double> edges_to_peri(N, 0.0);

      double selfloop = 0;
      for (auto& adj : G.neighbours(node_id)) {
        if (node_id == adj.node) {
          selfloop += adj.weight;
          continue;
        }

        edges_to_core[c[adj.node]] += adj.weight * x[adj.node];
        edges_to_peri[c[adj.node]] += adj.weight * (1 - x[adj.node]);
      }

      double D_core = sum_of_deg_core[c[node_id]] - deg * x[node_id];
      double D_peri = sum_of_deg_peri[c[node_id]] - deg * (1 - x[node_id]);
      double dQold = _calc_dQ_conf(edges_to_core[c[node_id]], edges_to_peri[c[node_id]], deg, D_core, D_peri,
                                   selfloop, x[node_id], resol);

      dQ = 0;
      for (auto& adj : G.neighbours(node_id)) {
        int cid = c[adj.node];

        D_core = sum_of_deg_core[cid] - deg * x[node_id] * (double)!!(c[node_id] == cid);
        D_peri = sum_of_deg_peri[cid] - deg * (1 - x[node_id]) * (double)!!(c[node_id] == cid);

        double Q_i_core =
            _calc_dQ_conf(edges_to_core[cid], edges_to_peri[cid], deg, D_core, D_peri, selfloop, 1, resol);
        double Q_i_peri =
            _calc_dQ_conf(edges_to_core[cid], edges_to_peri[cid], deg, D_core, D_peri, selfloop, 0, resol);
        Q_i_core -= dQold;
        Q_i_peri -= dQold;

        if (max(Q_i_core, Q_i_peri) < dQ) continue;

        if (Q_i_peri < Q_i_core) {
          xprime = 1;
          cprime = cid;
          dQ = Q_i_core;
        } else if (Q_i_peri > Q_i_core) {
          xprime = 0;
          cprime = cid;
          dQ = Q_i_peri;
        } else {
          cprime = cid;
          if (_udist(mtrnd) < 0.5) {
            xprime = 1;
          } else {
            xprime = 0;
          }
          dQ = Q_i_core;
        }
      }
    };

    int N = G.get_num_nodes();
    vector<double> sum_of_deg_core(N, 0.0);
    vector<double> sum_of_deg_peri(N, 0.0);
    vector<int> order(N);
    bool isupdated = false;
    c.clear();
    x.clear();
    c.assign(N, 0);
    x.assign(N, 1);
    for (int i = 0; i < N; i++) {
      order[i] = i;
      c[i] = i;
      sum_of_deg_core[i] += x[i] * phi[i];
    };

    do {
      isupdated = false;
      shuffle(order.begin(), order.end(), mtrnd);

      for (int scan_count = 0; scan_count < N; scan_count++) {
        int i = order[scan_count];

        int cprime = c[i];     // c'
        double xprime = x[i];  // x'

        double dQ = 0;
        double deg = phi[i];
        _propose_new_label(G, c, x, deg, sum_of_deg_core, sum_of_deg_peri, resol, i, cprime, xprime, dQ,
                           mtrnd);

        if (dQ <= 0) continue;

        if ((c[i] == cprime) & (x[i] == xprime)) continue;

        sum_of_deg_core[c[i]] -= deg * x[i];
        sum_of_deg_peri[c[i]] -= deg * (1 - x[i]);

        sum_of_deg_core[cprime] += deg * xprime;
        sum_of_deg_peri[cprime] += deg * (1 - xprime);

        c[i] = cprime;
        x[i] = xprime;

        isupdated = true;
      }

    } while (isupdated == true);

    /* Remove empty core-periphery pairs */
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

  /* Run the label switching algorithm num_run times and select
     the result yielding the largest quality
   */
  int N = G.get_num_nodes();
  c.clear();
  x.clear();
  c.assign(N, 0);
  x.assign(N, 1.0);

  Q = -1;
  for (int i = 0; i < num_of_runs; i++) {
    vector<int> ci;
    vector<double> xi;
    vector<double> qi;
    double Qi = 0.0;

    __label_switching(G, ci, xi, phi, resol, _mtrnd);

    _calc_Q_projected(G, ci, xi, phi, resol, qi);
    Qi = accumulate(qi.begin(), qi.end(), 0.0);

    if (Qi > Q) {
      c = ci;
      x = xi;
      q.clear();
      q = qi;
      Q = Qi;
    }
  }
}

/* Louvain algorithm */
void KM_multiresol::_louvain(Graph& G,
                             const int num_of_runs,
                             vector<int>& c,
                             vector<double>& x,
                             double& Q,
                             vector<double>& q,
                             map<int, int>& deg_port,
                             double resol,
                             mt19937_64& mtrnd) {
  /* Coarse-grain networks */
  auto _coarse_graining = [](Graph& G, const vector<int>& c, const vector<double>& x, Graph& newG,
                             vector<int>& toLayerId) {
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
    int N = c.size();
    vector<int> ids(N, 0);
    int maxid = 0;
    for (int i = 0; i < N; i++) {
      ids[i] = 2 * c[i] + (int)x[i];
      maxid = max(maxid, ids[i]);
    }
    _relabeling(ids);
    toLayerId.clear();
    toLayerId.assign(maxid + 1, 0);
    for (int i = 0; i < N; i++) {
      toLayerId[2 * c[i] + (int)x[i]] = ids[i];
    }

    newG = Graph();
    for (auto& node : G.adjacency_list()) {
      int i = node.first;
      int mi = 2 * c[i] + (int)x[i];
      for (auto& adj : node.second) {
        int mj = 2 * c[adj.node] + (int)x[adj.node];

        int sid = toLayerId[mi];
        int did = toLayerId[mj];
        newG.addEdge(sid, did, adj.weight);
      }
    }

    newG.aggregate_multi_edges();
  };

  /* Initialisation */
  int N = G.get_num_nodes();
  c.clear();
  x.clear();
  c.assign(N, 0);
  x.assign(N, 1);
  for (int i = 0; i < N; i++) c[i] = i;

  vector<int> ct = c;
  vector<double> xt = x;
  Graph cnet_G;
  vector<int> toLayerId;
  _coarse_graining(G, ct, xt, cnet_G, toLayerId);

  Q = 0;

  /* Main loop */
  int cnet_N;
  do {
    cnet_N = cnet_G.get_num_nodes();

    /* First step (Label switching) */
    vector<int> cnet_c;
    vector<double> cnet_x;
    double Qt = 0;
    vector<double> qt;
    _label_switching(cnet_G, num_of_runs, cnet_c, cnet_x, Qt, qt, deg_port, resol, mtrnd);

    for (int i = 0; i < N; i++) {
      int cnet_id = toLayerId[2 * ct[i] + xt[i]];
      ct[i] = cnet_c[cnet_id];
      xt[i] = cnet_x[cnet_id];
    }

    _calc_Q_projected(cnet_G, cnet_c, cnet_x, deg_port, resol, qt);
    Qt = accumulate(qt.begin(), qt.end(), 0.0);

    if (Qt >= Q) {
      c = ct;
      x = xt;
      Q = Qt;
      q = qt;
    }

    /* Second step (Coarse-graining) */
    Graph new_cnet_G;
    _coarse_graining(cnet_G, cnet_c, cnet_x, new_cnet_G, toLayerId);
    cnet_G = new_cnet_G;

    int sz = cnet_G.get_num_nodes();
    if (sz == cnet_N) break;

  } while (true);

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
}

void KM_multiresol::_est_null_model(Graph& G, vector<int>& ports, vector<int>& routes, int num_of_net) {
  auto _Chung_Lu_Algorithm_bipartite = [](const vector<double>& rdeg, const vector<double>& cdeg,
                                          const vector<int>& rnodes, const vector<int>& cnodes, Graph& G,
                                          bool isunweighted, mt19937_64& mtrnd) {

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
              G.addEdge(rnodes[u], cnodes[v] + Nr, w);
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
  for (auto& p : ports) rdeg[idx] = G.degree(p);
  idx = 0;
  for (auto& r : routes) cdeg[idx] = G.degree(r);

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
    Graph G_rand;
    vector<int> cr;
    vector<double> xr;
    double Qr;
    vector<double> qr;

    int tid = omp_get_thread_num();
    mt19937_64 mtrnd = mtrnd_list[tid];

    _Chung_Lu_Algorithm_bipartite(rdeg, cdeg, rdeg_rank, cdeg_rank, G_rand, isunweighted, mtrnd);

    Graph G_rand_projected;
    _one_mode_projection(G_rand, ports, routes, _phi, G_rand_projected);

    _louvain(G_rand_projected, _num_runs, cr, xr, Qr, qr, _phi, _resol, mtrnd);

    int K_rand = qr.size();
    vector<int> nsr(K_rand, 0);
    for (int i = 0; i < Nport; i++) {
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

vector<double> KM_multiresol::_compute_p_values(Graph& G,
                                                vector<int>& c,
                                                vector<double>& x,
                                                vector<int>& ports,
                                                vector<int>& routes,
                                                vector<double>& phi,
                                                double resol) {
  auto normcdf = [](double x) { return 0.5 + 0.5 * erf(x * M_SQRT1_2); };

  int K = *max_element(c.begin(), c.end()) + 1;

  vector<double> q;
  calc_Q(G, c, x, ports, routes, phi, resol, q);

  vector<double> n(K, 0.0);
  for (auto& ci : c) n[ci]++;

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
      double qbar = _qhat[s] + sig_nq / sig_nn * (double)(n[k] - _nhat[s]);

      double t = sig_nn * (q[k] - qbar) / (sqrt(sig_nn * sig_qq - sig_nq * sig_nq) * h);
      double cum = normcdf(t);

      double w = exp(-(double)pow(n[k] - _nhat[s], 2) / (2.0 * h * h * sig_nn)) + 1e-33;
      numer += cum * w;
      denom += w;
    }
    p_values[k] = 1.0 - numer / denom;
  }
  return p_values;
}

void KM_multiresol::_one_mode_projection(
    Graph& G, vector<int>& ports, vector<int>& routes, vector<double>& phi, Graph& G_projected) {
  for (auto& r : routes) {
    vector<AdjacentNode> adj = G.neighbours(r);
    int sz = adj.size();
    for (int i = 0; i < sz; i++) {
      for (int j = i + 1; j < sz; j++) {
        G_projected.addEdge(adj[i].node, adj[j].node, phi[r] / (double)(G.degree(r) - 1));
      }
    }
  }

  G_projected.aggregate_multi_edges();
}

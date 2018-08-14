#include <omp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include <graph.h>

using namespace std;

class KM_multiresol {
 public:
  // Constructor
  KM_multiresol();
  KM_multiresol(int num_runs);

  // detect and calc_Q take bipartite networks
  void detect(Graph& G, vector<double>& theta, double resol);

  void calc_Q(Graph& G,  // bipartite network
              const vector<int>& c,
              const vector<double>& x,
              vector<double>& theta,
              double resol,
              vector<double>& q);

  // Getter
  vector<int> get_c() const { return _c; };
  vector<double> get_x() const { return _x; };
  vector<double> get_q() const { return _q; };

 private:
  /* For core-periphery detection */
  int _num_runs;
  vector<int> _c;
  vector<double> _x;
  double _Q;
  vector<double> _q;
  vector<double> _theta;
  double _resol;

  /* Random number generator */
  mt19937_64 _mtrnd;

  void _label_switching(Graph& G,
                        const int num_of_runs,
                        vector<int>& c,
                        vector<double>& x,
                        double& Q,
                        vector<double>& q,
                        vector<double>& theta,
                        double resol,
                        mt19937_64& mtrnd);

  void _louvain(Graph& G,
                const int num_of_runs,
                vector<int>& c,
                vector<double>& x,
                double& Q,
                vector<double>& q,
                vector<double>& theta,
                double resol,
                mt19937_64& mtrnd);
};

/*-----------------------------
Constructor
-----------------------------*/

KM_multiresol::KM_multiresol(int num_runs) {
  mt19937_64 mtrnd;
  random_device r;
  seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
  mtrnd.seed(seed);
  _mtrnd = mtrnd;

  _num_runs = num_runs;
};

KM_multiresol::KM_multiresol() { KM_multiresol(10); };

/*-----------------------------
Public functions
-----------------------------*/
void KM_multiresol::detect(Graph& G, vector<double>& theta, double resol) {
  _theta = theta;
  _resol = resol;
  _louvain(G, _num_runs, _c, _x, _Q, _q, theta, resol, _mtrnd);
  //_label_switching(G, _num_runs, _c, _x, _Q, _q, theta, resol, _mtrnd);
}

void KM_multiresol::calc_Q(Graph& G,
                           const vector<int>& c,
                           const vector<double>& x,
                           vector<double>& theta,
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
    Dc[c[i]] += x[i] * theta[i];
    Dp[c[i]] += (1 - x[i]) * theta[i];
  }
  Omega /= 2;

  for (int k = 0; k < K; k++) {
    q[k] = (q[k] - resol * (Dc[k] * Dc[k] + 2 * Dc[k] * Dp[k])) / (2 * Omega);
  }
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
                                     vector<double>& theta,
                                     double resol,
                                     mt19937_64& mtrnd) {
  /* Label switching algorithm */
  auto __label_switching = [](Graph& G, vector<int>& c, vector<double>& x, vector<double>& theta,
                              double resol, mt19937_64& mtrnd) {

    auto _propose_new_label = [](Graph& G, const vector<int>& c, const vector<double>& x, const double th,
                                 const vector<double>& sum_of_th_core, const vector<double>& sum_of_th_peri,
                                 const double resol, const int node_id, int& cprime, double& xprime,
                                 double& dQ, mt19937_64& mtrnd) {

      auto _calc_dQ_conf = [](double d_i_c, double d_i_p, double d_i, double D_c, double D_p, double selfloop,
                              double x, const double resol) {
        return 2 * (d_i_c + d_i_p * (x)-resol * d_i * (D_c + D_p * x)) +
               x * (selfloop - 2 * resol * d_i * d_i);
      };  // End of calc_dQ_conf

      cprime = c[node_id];
      xprime = x[node_id];

      if (G.isnode(node_id) == false) return;

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

      double D_core = sum_of_th_core[c[node_id]] - th * x[node_id];
      double D_peri = sum_of_th_peri[c[node_id]] - th * (1 - x[node_id]);
      double dQold = _calc_dQ_conf(edges_to_core[c[node_id]], edges_to_peri[c[node_id]], th, D_core, D_peri,
                                   selfloop, x[node_id], resol);
      dQ = 0;
      for (auto& adj : G.neighbours(node_id)) {
        int cid = c[adj.node];

        D_core = sum_of_th_core[cid] - th * x[node_id] * (double)!!(c[node_id] == cid);
        D_peri = sum_of_th_peri[cid] - th * (1 - x[node_id]) * (double)!!(c[node_id] == cid);

        double Q_i_core =
            _calc_dQ_conf(edges_to_core[cid], edges_to_peri[cid], th, D_core, D_peri, selfloop, 1, resol);
        double Q_i_peri =
            _calc_dQ_conf(edges_to_core[cid], edges_to_peri[cid], th, D_core, D_peri, selfloop, 0, resol);
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
    };  // End of _propose_new_label

    int N = G.get_num_nodes();
    vector<double> sum_of_th_core(N, 0.0);
    vector<double> sum_of_th_peri(N, 0.0);
    vector<int> order(N);
    bool isupdated = false;
    c.clear();
    x.clear();
    c.assign(N, 0);
    x.assign(N, 1);
    for (int i = 0; i < N; i++) {
      order[i] = i;
      c[i] = i;
      sum_of_th_core[i] += x[i] * theta[i];
    };

    do {
      isupdated = false;
      shuffle(order.begin(), order.end(), mtrnd);
      for (int scan_count = 0; scan_count < N; scan_count++) {
        int i = order[scan_count];

        int cprime = c[i];     // c'
        double xprime = x[i];  // x'

        double dQ = 0;
        double th = theta[i];
        _propose_new_label(G, c, x, th, sum_of_th_core, sum_of_th_peri, resol, i, cprime, xprime, dQ, mtrnd);

        if (dQ <= 0) continue;

        if ((c[i] == cprime) & (x[i] == xprime)) continue;

        sum_of_th_core[c[i]] -= th * x[i];
        sum_of_th_peri[c[i]] -= th * (1 - x[i]);

        sum_of_th_core[cprime] += th * xprime;
        sum_of_th_peri[cprime] += th * (1 - xprime);

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
  };  // End of __label_switching

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

    __label_switching(G, ci, xi, theta, resol, mtrnd);

    vector<double> qi;
    calc_Q(G, ci, xi, theta, resol, qi);
    double Qi = accumulate(qi.begin(), qi.end(), 0.0);
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
                             vector<double>& theta,
                             double resol,
                             mt19937_64& mtrnd) {
  /* Coarse-grain networks */
  auto _coarse_graining = [](Graph& G, const vector<int>& c, const vector<double>& x, Graph& newG,
                             vector<int>& toLayerId, vector<double>& theta, vector<double>& theta_new) {
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
      int sid = toLayerId[mi];
      for (auto& adj : node.second) {
        if (i > adj.node) continue;
        int mj = 2 * c[adj.node] + (int)x[adj.node];
        int did = toLayerId[mj];

        if (sid == did) {
          newG.addEdge(sid, did, adj.weight / 2);
        } else {
          newG.addEdge(sid, did, adj.weight);
        }
      }
    }
    newG.aggregate_multi_edges();

    theta_new.clear();
    theta_new.assign(newG.get_num_nodes(), 0);
    for (auto& node : G.adjacency_list()) {
      int i = node.first;
      int mi = 2 * c[i] + (int)x[i];
      int sid = toLayerId[mi];
      theta_new[sid] += theta[i];
    }
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
  vector<double> cnet_theta;
  _coarse_graining(G, ct, xt, cnet_G, toLayerId, theta, cnet_theta);

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
    // cout<<cnet_G.get_total_edge_weight()<<" "<<G.get_total_edge_weight()<<endl;
    _label_switching(cnet_G, num_of_runs, cnet_c, cnet_x, Qt, qt, cnet_theta, resol, mtrnd);

    for (int i = 0; i < N; i++) {
      int cnet_id = toLayerId[2 * ct[i] + xt[i]];
      ct[i] = cnet_c[cnet_id];
      xt[i] = cnet_x[cnet_id];
    }

    calc_Q(cnet_G, cnet_c, cnet_x, cnet_theta, resol, qt);
    // calc_Q(G, ct, xt, theta, resol, qt);
    Qt = accumulate(qt.begin(), qt.end(), 0.0);

    if (Qt >= Q) {
      c = ct;
      x = xt;
      Q = Qt;
      q = qt;
    }

    /* Second step (Coarse-graining) */
    Graph new_cnet_G;
    vector<double> new_cnet_theta;
    _coarse_graining(cnet_G, cnet_c, cnet_x, new_cnet_G, toLayerId, cnet_theta, new_cnet_theta);
    cnet_G = new_cnet_G;
    cnet_theta = new_cnet_theta;

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

#include <omp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <vector>

#include <km_multiresol.h>

using namespace std;

class KM_omp {
 public:
  // Constructor
  KM_omp();

  void detect(Graph& biG, vector<int>& ports, vector<int>& routes, map<int, double>& phi, double resol);

  /* Getters */
  vector<int> get_c() const { return _c; };
  vector<double> get_x() const { return _x; };

  /* Setters */
  void set_num_of_runs(int param) { _num_runs_KM_multiresol = param; };
  void set_significance_level(int param) { _significance_level = param; };
  void set_num_of_samples(int param) { _num_samples = param; };
  void set_num_of_rand_nets(int param) { _num_rand_nets = param; };
  void set_consensus_threshold(double param) { _consensus_th = param; };

 private:
  int _num_runs_KM_multiresol;
  int _num_samples;
  int _num_rand_nets;
  double _significance_level;
  double _consensus_th;

  vector<double> _nhat;
  vector<double> _qhat;
  vector<int> _c;
  vector<double> _x;

  void _one_mode_projection(
      Graph& biG, vector<int>& ports, vector<int>& routes, map<int, double>& phi, Graph& uniG);

  void _est_null_model(Graph& biG,
                       vector<int>& ports,
                       vector<int>& routes,
                       map<int, double>& phi,
                       double resol,
                       int num_of_net);

  vector<double> _calc_p_values(vector<double>& q, vector<double>& n);
  vector<int> _connected_components(Graph& U, double th, int N);
};

KM_omp::KM_omp() {
  _num_runs_KM_multiresol = 10;
  _num_samples = 100;
  _num_rand_nets = 500;
  _significance_level = 0.05;
  _consensus_th = 0.9;
}

void KM_omp::detect(
    Graph& biG, vector<int>& ports, vector<int>& routes, map<int, double>& phi, double resol) {
  int Np = ports.size();
  int Nr = routes.size();

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

  Graph uniG;
  _one_mode_projection(biG, ports, routes, phi, uniG);

  _est_null_model(biG, ports, routes, phi, resol, _num_rand_nets);

  Graph U;
  vector<double> X(Np, 0.0);
  vector<double> sig_count(Np, 0.0);
#ifdef _OPENMP
#pragma omp parallel for shared(X, sig_count, U)
#endif
  for (int sid = 0; sid < _num_samples; sid++) {
    KM_multiresol _km(_num_runs_KM_multiresol);
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
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      for (int i = 0; i < Np; i++) {
        if (pvals[cs[i]] > alpha) continue;

        for (int j = i + 1; j < Np; j++) {
          if (pvals[cs[j]] > alpha) continue;

          if (cs[i] == cs[j]) {
            U.addEdge(i, j, 1);
          }
        }
        X[i] += xs[i];
        sig_count[i]++;
      }
    }
  }
  for (int i = 0; i < Np; i++) {
    if (sig_count[i] == 0) continue;
    X[i] /= sig_count[i];
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
    double phi_r = phi[r];
    double dr = G.degree(r);
    double w = phi_r / (double)(dr - 1.0);
    for (int i = 0; i < sz; i++) {
      for (int j = i + 1; j < sz; j++) {
        uniG.addEdge(adj[i].node, adj[j].node, w);
      }
    }
  }
  uniG.aggregate_multi_edges();
}

void KM_omp::_est_null_model(Graph& biG,
                             vector<int>& ports,
                             vector<int>& routes,
                             map<int, double>& phi,
                             double resol,
                             int num_of_net) {
  auto _Chung_Lu_Algorithm_bipartite = [](const vector<double>& rdeg, const vector<double>& cdeg,
                                          const vector<int>& rnodes, const vector<int>& cnodes,
                                          vector<int>& rnodes_ids, vector<int>& cnodes_ids, Graph& biG,
                                          bool isunweighted, mt19937_64& mtrnd) {

    uniform_real_distribution<double> udist(0.0, 1.0);
    int Nr = rdeg.size();
    int Nc = cdeg.size();
    double M = accumulate(cdeg.begin(), cdeg.end(), 0.0);

    vector<vector<pair<int, double>>> tmp(Nr, vector<pair<int, double>>(0));
    for (int u = 0; u < Nr; u++) {
      for (int v = 0; v < Nc; v++) {
        double p = min(1.0, rdeg[rnodes[u]] * cdeg[cnodes[v]] / M);
        while (v < Nc && p > 0) {
          if (p != 1) {
            double r = udist(mtrnd);
            v = v + floor(log(r) / log(1 - p));
          }
          if (v < Nc) {
            double q = min(rdeg[rnodes[u]] * cdeg[cnodes[v]] / M, 1.0);
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
              biG.addEdge(rnodes_ids[rnodes[u]], cnodes_ids[cnodes[v]], 1);
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

  vector<int> rdeg_rank(Nport);
  iota(rdeg_rank.begin(), rdeg_rank.end(), 0);
  sort(rdeg_rank.begin(), rdeg_rank.end(), [&](int x, int y) { return rdeg[x] > rdeg[y]; });
  vector<int> cdeg_rank(Nroute);
  iota(cdeg_rank.begin(), cdeg_rank.end(), 0);
  sort(cdeg_rank.begin(), cdeg_rank.end(), [&](int x, int y) { return cdeg[x] > cdeg[y]; });

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
#pragma omp parallel for shared(mtrnd_list, nhat, qhat)
#endif
  for (int it = 0; it < _num_rand_nets; it++) {
    int tid = omp_get_thread_num();

    Graph biG_rand;
    Graph uniG_rand;
    _Chung_Lu_Algorithm_bipartite(rdeg, cdeg, rdeg_rank, cdeg_rank, ports, routes, biG_rand, isunweighted,
                                  mtrnd_list[tid]);

    _one_mode_projection(biG_rand, ports, routes, phi, uniG_rand);
    map<int, int> node2node = uniG_rand.renumbering();
    vector<double> theta(uniG_rand.get_num_nodes(), 0.0);
    for (auto& p : uniG_rand.adjacency_list()) {
      theta[p.first] = biG_rand.degree(ports[node2node[p.first]]);
    }

    KM_multiresol _km(_num_runs_KM_multiresol);
    _km.detect(uniG_rand, theta, resol);
    vector<int> cr = _km.get_c();
    vector<double> qr = _km.get_q();

    int K_rand = qr.size();
    vector<int> nsr(K_rand, 0);
    for (auto& p : uniG_rand.adjacency_list()) {
      nsr[cr[p.first]]++;
    }
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      for (int k = 0; k < K_rand; k++) {
        if (nsr[k] <= 1) continue;
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
  double q_ave = (double)accumulate(_nhat.begin(), _nhat.end(), 0.0) / (double)S;
  double s_ave = (double)accumulate(_qhat.begin(), _qhat.end(), 0.0) / (double)S;
  double s_std = 0;
  double q_std = 0;
  for (int s = 0; s < S; s++) {
    s_std += pow((double)_nhat[s] - s_ave, 2);
    q_std += pow(_qhat[s] - q_ave, 2);
  }
  s_std = sqrt(s_std / (double)(S - 1));
  q_std = sqrt(q_std / (double)(S - 1));

  double gamma = 0;
  if ((s_std <= 1e-30) | (q_std <= 1e-30)) {
    s_std = 1e-20;
  } else {
    for (int s = 0; s < S; s++) {
      gamma += ((double)_nhat[s] - s_ave) * (_qhat[s] - q_ave);
    }
    gamma = gamma / ((double)(S - 1) * s_std * q_std);
  }

  double h = max(pow((double)S, -1.0 / 6.0), 1e-32);
  vector<double> p_values(K, 1.0);
  for (int k = 0; k < K; k++) {
    if ((s_std <= 1e-30) | (q_std <= 1e-30)) {
      continue;
    }
    double numer = 0.0;
    double denom = 0.0;
    for (int s = 0; s < S; s++) {
      double w = exp(-(double)pow(n[k] - _nhat[s], 2) / (2.0 * h * h * s_std * s_std));
      double cum = normcdf(((q[k] - _qhat[s]) / (h * q_std) - gamma * (n[k] - _nhat[s]) / (h * s_std)) /
                           sqrt(1.0 - gamma * gamma));

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
          xnew[nodes.first] += x[nodes.first] + x[adj.node];
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
  return c;
}

#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <iostream>

#include <km_multiresol.h>
#include <km_omp.h>

using namespace std;
namespace py = pybind11;

py::list _detect(py::array_t<int> edges,
                 py::array_t<int> ports,
                 py::array_t<int> routes,
                 py::array_t<double> phi,
                 double resol,
                 int num_samples,
                 int num_runs,
                 double consensus_threshold,
                 double significance_level,
                 int num_rand_nets) {

  auto edges2graph = [](py::array_t<int> edges_array_t) {
    Graph G;
    auto edges = edges_array_t.data();
    auto r = edges_array_t.request();
    int M = r.shape[0];
    for (int i = 0; i < M; i++) {
      int sid = edges[2 * i];
      int did = edges[2 * i + 1];
      G.addEdge(sid, did, 1);
    }
    G.aggregate_multi_edges();
    return G;
  };

  auto pyarray2vec_int = [](py::array_t<int> a) {
    auto r = a.request();
    int N = r.shape[0];
    vector<int> a_vec(N);
    auto a_py = a.data();
    for (int i = 0; i < N; i++) a_vec[i] = a_py[i];
    return a_vec;
  };

  auto pyarray2phi_map = [](py::array_t<int> rt, py::array_t<double> ph) {
    auto rt_py = rt.data();
    auto ph_py = ph.data();
    auto r = rt.request();
    int sz = r.shape[0];
    map<int, double> retval;
    for (int i = 0; i < sz; i++) {
      retval.insert(make_pair(rt_py[i], ph_py[i]));
    }
    return retval;
  };

  auto cx2list = [](vector<int>& c, vector<double>& x) {
    py::list results(2);
    int N = c.size();
    py::array_t<double> cids_array_t(N);
    auto cids = cids_array_t.mutable_data();
    py::array_t<double> xs_array_t(N);
    auto xs = xs_array_t.mutable_data();
    for (int i = 0; i < N; i++) {
      cids[i] = c[i];
      xs[i] = x[i];
    }
    results[0] = cids_array_t;
    results[1] = xs_array_t;
    return results;
  };

  /* Construct a graph object */
  Graph G = edges2graph(edges);

  /* Convert python objects to c++ objects */
  map<int, double> phi_map = pyarray2phi_map(routes, phi);
  vector<int> ports_vec = pyarray2vec_int(ports);
  vector<int> routes_vec = pyarray2vec_int(routes);
  /* Core-periphery detection */
  KM_omp km = KM_omp();
  km.set_num_of_runs(num_runs);
  km.set_significance_level(significance_level);
  km.set_num_of_samples(num_samples);
  km.set_num_of_rand_nets(num_rand_nets);
  km.set_consensus_threshold(consensus_threshold);
  km.detect(G, ports_vec, routes_vec, phi_map, resol);

  /* Retrieve results */
  vector<int> c = km.get_c();
  vector<double> x = km.get_x();

  /* Convert c++ objects to python objects */
  py::list results = cx2list(c, x);
  return results;
}

PYBIND11_MODULE(_km_ompnet, m) {
  m.doc() = "The KM algorithm for networks constructed from one-mode-projection. ";

  // CP detection algorithm
  m.def("_detect", &_detect, "KM algorithm", py::arg("edges"), py::arg("ports"), py::arg("routes"),
        py::arg("phi"), py::arg("resol"), py::arg("num_samples"), py::arg("num_runs"),
        py::arg("consensus_threshold"), py::arg("significance_level"),
        py::arg("num_rand_nets"));
}

#include <bilad.h>
#include <cassert>
#include <cmath>
#include <limits>

using namespace BiLAD;
using namespace std;

namespace BiLAD {
tuple<GraphPath, GraphPath, size_t> PSQSR(const Graph& graph, size_t src,
                                          size_t dest, weight_t delta) {
    vector<double> c_dist(graph.size(), INFINITY);
    vector<double> d_dist(graph.size(), INFINITY);
    vector<size_t> prev(graph.size(), -1);
    size_t count = 0;

    c_dist[src] = 0;
    d_dist[src] = 0;

    for (size_t i = 0; i < graph.size(); ++i) {
        for (size_t j = 0; j < graph.size(); ++j) {
            for (size_t k = 0; k < graph[j].size(); ++k) {
                const auto& edge = graph[j][k];
                size_t u = j;
                size_t v = edge.first;
                auto weight = edge.second;

                auto x_u = c_dist[u];
                auto y_u = d_dist[u];
                auto x_v = c_dist[v];
                auto y_v = d_dist[v];
                if (isinf(x_u) || isinf(y_u))
                    continue;
                if (isinf(x_v) || isinf(y_v)) {
                    c_dist[v] = x_u + weight.first;
                    d_dist[v] = y_u + weight.second;
                    prev[v] = u;
                    continue;
                }

                auto p = x_u + weight.first - x_v;
                auto q = y_u + weight.second - y_v;
                if (p >= 0 && q >= 0) {
                    continue;
                } else if (p <= 0 && q <= 0) {
                    c_dist[v] = x_u + weight.first;
                    d_dist[v] = y_u + weight.second;
                    prev[v] = u;
                    continue;
                }
                auto lambda = -p / q;
                count++;
                auto [p_c, p_d] = biweight_dijkstra(graph, src, dest, lambda);

                if (graph.d_dist(p_d) > delta) {
                    if (q < 0) {
                        c_dist[v] = x_u + weight.first;
                        d_dist[v] = y_u + weight.second;
                        prev[v] = u;
                        continue;
                    } else {
                        continue;
                    }
                } else if (graph.d_dist(p_c) > delta) {
                    if (q < 0) {
                        continue;
                    } else {
                        c_dist[v] = x_u + weight.first;
                        d_dist[v] = y_u + weight.second;
                        prev[v] = u;
                        continue;
                    }
                } else if (graph.d_dist(p_c) > delta &&
                           graph.d_dist(p_d) < delta) {
                    return make_tuple(p_d, p_c, count);
                } else {
                    assert(false);
                }
            }
        }
    }

    GraphPath path;
    size_t i = dest;
    path.push_back(dest);
    while (prev[i] != src) {
        path.push_back(prev[i]);
        i = prev[i];
    }
    path.push_back(src);
    reverse(path.begin(), path.end());

    return make_tuple(path, GraphPath{}, count);
}
} // namespace BiLAD
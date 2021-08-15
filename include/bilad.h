#ifndef BILAD_BILAD_H
#define BILAD_BILAD_H

#include <algorithm>
#include <cstdint>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

// TODO: Change API to multiple Adj Table

namespace BiLAD {
const double EPSILON = 1e-5;
using weight_t = double;
using GraphPath = std::vector<size_t>;

class Graph {
  private:
    using AdjTable = std::vector<
        std::vector<std::pair<size_t, std::pair<weight_t, weight_t>>>>;
    AdjTable adj;

  public:
    explicit Graph(size_t n) : adj(n) {}

    void add_edge(size_t i, size_t j, weight_t c, weight_t d) {
        adj[i].emplace_back(std::make_pair(j, std::make_pair(c, d)));
        adj[j].emplace_back(std::make_pair(i, std::make_pair(c, d)));
    }

    void add_directed_edge(size_t i, size_t j, weight_t c, weight_t d) {
        adj[i].emplace_back(std::make_pair(j, std::make_pair(c, d)));
    }

    void delete_edge(size_t i, size_t j) {
        {
            auto pos1 =
                std::find_if(adj[i].begin(), adj[i].end(),
                             [&](const auto& v) { return v.first == j; });
            if (pos1 != adj[i].end()) {
                std::swap(*pos1, adj[i].back());
                adj[i].pop_back();
            }
        }

        {
            auto pos2 =
                std::find_if(adj[j].begin(), adj[j].end(),
                             [&](const auto& v) { return v.first == i; });
            if (pos2 != adj[j].end()) {
                std::swap(*pos2, adj[j].back());
                adj[j].pop_back();
            }
        }
    }

    void delete_vertex(const std::unordered_set<size_t>& d) {
        for (auto& edge : adj) {
            for (size_t i = 0; i < edge.size(); ++i) {
                for (auto v : d) {
                    auto pos = std::find_if(
                        edge.begin(), edge.end(),
                        [&](const auto& val) { return val.first == v; });
                    if (pos != edge.end())
                        edge.erase(pos);
                }
            }
        }
        for (auto v : d) {
            adj[v].clear();
        }
    }

    auto& operator[](size_t i) const { return adj[i]; }

    size_t size() const { return adj.size(); }

    weight_t d_dist(const GraphPath& path) const {
        weight_t dist = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            const auto& edges = adj[path[i - 1]];
            const auto& target =
                *std::find_if(edges.begin(), edges.end(),
                              [&](auto v) { return v.first == path[i]; });
            dist += target.second.second;
        }
        return dist;
    }

    weight_t c_dist(const GraphPath& path) const {
        weight_t dist = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            const auto& edges = adj[path[i - 1]];
            const auto& target =
                *std::find_if(edges.begin(), edges.end(),
                              [&](auto v) { return v.first == path[i]; });
            dist += target.second.first;
        }
        return dist;
    }

    double weighted_dist(const GraphPath& path, double lambda,
                         double gama = 1) const {
        double dist = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            const auto& edges = adj[path[i - 1]];
            const auto& target =
                *std::find_if(edges.begin(), edges.end(),
                              [&](auto v) { return v.first == path[i]; });
            dist += gama * target.second.first + lambda * target.second.second;
        }
        return dist;
    }

    void preprocess() {
        for (auto& Edge : adj) {
            for (auto& edge : Edge) {
                edge.second.first *= 2;
                edge.second.second *= 2;
            }
        }
    }
};

class TriWeightedGraph {
  private:
    using AdjTable = std::vector<std::vector<
        std::pair<size_t, std::tuple<weight_t, weight_t, weight_t>>>>;
    AdjTable adj;

  public:
    explicit TriWeightedGraph(size_t n) : adj(n) {}

    void add_edge(size_t i, size_t j, weight_t c, weight_t d1, weight_t d2) {
        adj[i].emplace_back(std::make_pair(j, std::make_tuple(c, d1, d2)));
        adj[j].emplace_back(std::make_pair(i, std::make_tuple(c, d1, d2)));
    }

    void add_directed_edge(size_t i, size_t j, weight_t c, weight_t d1,
                           weight_t d2) {
        adj[i].emplace_back(std::make_pair(j, std::make_tuple(c, d1, d2)));
    }

    auto& operator[](size_t i) const { return adj[i]; }

    size_t size() const { return adj.size(); }

    weight_t d2_dist(const GraphPath& path) const {
        weight_t dist = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            const auto& edges = adj[path[i - 1]];
            const auto& target =
                *std::find_if(edges.begin(), edges.end(),
                              [&](auto v) { return v.first == path[i]; });
            dist += std::get<2>(target.second);
        }
        return dist;
    }

    weight_t d1_dist(const GraphPath& path) const {
        weight_t dist = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            const auto& edges = adj[path[i - 1]];
            const auto& target =
                *std::find_if(edges.begin(), edges.end(),
                              [&](auto v) { return v.first == path[i]; });
            dist += std::get<1>(target.second);
        }
        return dist;
    }

    weight_t c_dist(const GraphPath& path) const {
        weight_t dist = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            const auto& edges = adj[path[i - 1]];
            const auto& target =
                *std::find_if(edges.begin(), edges.end(),
                              [&](auto v) { return v.first == path[i]; });
            dist += std::get<0>(target.second);
        }
        return dist;
    }

    void preprocess() {
        for (auto& Edge : adj) {
            for (auto& edge : Edge) {
                auto& [c, d1, d2] = edge.second;
                c *= 2;
                d1 *= 2;
                d2 *= 2;
            }
        }
    }

    Graph copy_c_d1() const {
        Graph g(this->size());
        for (size_t i = 0; i < adj.size(); ++i) {
            for (auto& edge : adj[i]) {
                auto& [c, d1, d2] = edge.second;
                g.add_directed_edge(i, edge.first, c, d1);
            }
        }
        return g;
    }

    Graph copy_c_d2() const {
        Graph g(this->size());
        for (size_t i = 0; i < adj.size(); ++i) {
            for (auto& edge : adj[i]) {
                auto& [c, d1, d2] = edge.second;
                g.add_directed_edge(i, edge.first, c, d2);
            }
        }
        return g;
    }

    Graph copy_d1_d2() const {
        Graph g(this->size());
        for (size_t i = 0; i < adj.size(); ++i) {
            for (auto& edge : adj[i]) {
                auto& [c, d1, d2] = edge.second;
                g.add_directed_edge(i, edge.first, d1, d2);
            }
        }
        return g;
    }

    void swap_d1_d2() {
        for (size_t i = 0; i < adj.size(); ++i) {
            for (auto& edge : adj[i]) {
                auto& [c, d1, d2] = edge.second;
                std::swap(d1, d2);
            }
        }
    }

    Graph copy_c_lambda_d1(double lambda) {
        Graph g(this->size());
        for (size_t i = 0; i < adj.size(); ++i) {
            for (auto& edge : adj[i]) {
                auto& [c, d1, d2] = edge.second;
                g.add_directed_edge(i, edge.first, c + lambda * d1, d2);
            }
        }
        return g;
    }
};

std::tuple<GraphPath, GraphPath> biweight_dijkstra(const Graph& graph,
                                                   size_t src, size_t dest,
                                                   double lambda,
                                                   double gama = 1);

std::tuple<GraphPath, GraphPath, double, bool, size_t>
bilad(const Graph& graph, size_t src, size_t dest, weight_t delta);

std::vector<GraphPath> yen_algorithm(const Graph& graph, size_t src,
                                     size_t dest, double lambda, size_t K,
                                     size_t* count_ptr = nullptr);

std::tuple<GraphPath, size_t> exact_bilad(const Graph& graph, size_t src,
                                          size_t dest, weight_t delta);

std::tuple<GraphPath, GraphPath, size_t> PSQSR(const Graph& graph, size_t src,
                                               size_t dest, weight_t delta);

typedef std::tuple<GraphPath, GraphPath, double, size_t>
SCSP(const Graph& graph, size_t src, size_t dest, weight_t delta);

std::tuple<GraphPath, GraphPath, GraphPath, GraphPath, double, double, size_t,
           size_t>
new_bilad(TriWeightedGraph& graph, size_t src, size_t dest, weight_t &delta1,
          weight_t &delta2, SCSP scsp);

static inline std::tuple<GraphPath, GraphPath, double, size_t>
bind_bilad(const Graph& graph, size_t src, size_t dest, weight_t delta) {
    auto [p_c, p_d, lambda, flag, count] = bilad(graph, src, dest, delta);
    return {p_c, p_d, lambda, count};
}

} // namespace BiLAD
#endif // BILAD_BILAD_H

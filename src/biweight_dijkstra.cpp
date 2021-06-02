#include <bilad.h>
#include <vector>
#include <limits>
#include <queue>
#include <utility>
#include <cmath>
#include <algorithm>

using namespace std;

namespace BiLAD {
    std::tuple<GraphPath, GraphPath>
    biweight_dijkstra(const Graph &graph, size_t src, size_t dest, double lambda, double gama) {
        size_t n = graph.size();

        vector<bool> visited(n, false);
        vector<double> dist(n, numeric_limits<double>::max());
        dist[src] = 0;
        vector<weight_t> c_dist(n, numeric_limits<weight_t>::max()), d_dist(n, numeric_limits<weight_t>::max());
        c_dist[src] = 0;
        d_dist[src] = 0;
        using Dist = pair<size_t, double>;
        auto cmp = [](const Dist &lhs, const Dist &rhs) { return lhs.second > rhs.second; };
        priority_queue<Dist, vector<Dist>, decltype(cmp)> pq(cmp);
        pq.push(make_pair(src, dist[src]));
        vector<size_t> c_prev(n), d_prev(n);

        while (!pq.empty()) {
            size_t curr = pq.top().first;
            pq.pop();
            if (visited[curr])continue;
            visited[curr] = true;
            for (const auto &i:graph[curr]) {
                size_t v = i.first;
                const auto &w = i.second;
                double weight = gama * w.first + lambda * w.second;
                double new_dist = dist[curr] + weight;
                if (dist[v] > new_dist) {
                    dist[v] = new_dist;
                    pq.push(make_pair(v, new_dist));

                    c_dist[v] = c_dist[curr] + w.first;
                    d_dist[v] = d_dist[curr] + w.second;

                    c_prev[v] = curr;
                    d_prev[v] = curr;
                } else if (abs(dist[v] - new_dist) < EPSILON) {
                    weight_t new_c_dist = c_dist[curr] + w.first;
                    if (c_dist[v] > new_c_dist) {
                        c_dist[v] = new_c_dist;
                        c_prev[v] = curr;
                    }

                    weight_t new_d_dist = d_dist[curr] + w.second;
                    if (d_dist[v] > new_d_dist) {
                        d_dist[v] = new_d_dist;
                        d_prev[v] = curr;
                    }
                }
            }
        }

        if (!visited[dest])
            return {{},
                    {}};

        GraphPath c_path, d_path;
        size_t i = dest;
        c_path.push_back(dest);
        while (c_prev[i] != src) {
            c_path.push_back(c_prev[i]);
            i = c_prev[i];
        }
        c_path.push_back(src);
        reverse(c_path.begin(), c_path.end());

        size_t j = dest;
        d_path.push_back(dest);
        while (d_prev[j] != src) {
            d_path.push_back(d_prev[j]);
            j = d_prev[j];
        }
        d_path.push_back(src);
        reverse(d_path.begin(), d_path.end());

        return {c_path, d_path};
    }
}
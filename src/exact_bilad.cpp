#include <bilad.h>
#include <queue>
#include <unordered_set>
#include <vector>

using namespace std;

namespace BiLAD {
tuple<GraphPath, size_t> exact_bilad(const Graph& graph, size_t src,
                                     size_t dest, weight_t delta) {
    auto [p1, p2, lambda, flag, count] = bilad(graph, src, dest, delta);
    if (!flag) {
        if (p1.empty() && p2.empty()) {
            return {GraphPath{}, count};
        } else if (p1.empty()) {
            return {p2, count};
        } else if (p2.empty()) {
            return {p1, count};
        } else {
            return {p2, count};
        }
    }
    auto p_star = p2;
    double L_B = graph.c_dist(p2) +
                 lambda * (graph.d_dist(p2) - static_cast<double>(delta));
    double U_B = graph.c_dist(p2);
    double gap = U_B - L_B;
    double tolerance = 1e-5;
    size_t k = 2;
    while (gap > tolerance) {
        k++;
        auto P = yen_algorithm(graph, src, dest, lambda, k, &count);
        if (P.size() != k)
            break;
        const GraphPath p_k = move(P.back());
        L_B = graph.c_dist(p_k) +
              lambda * (graph.d_dist(p_k) - static_cast<double>(delta));
        if (graph.d_dist(p_k) <= delta && graph.c_dist(p_k) < U_B) {
            p_star = p_k;
            U_B = graph.c_dist(p_star);
        }
        gap = U_B - L_B;
    }
    return {p_star, count};
}

vector<GraphPath> yen_algorithm(const Graph& graph, size_t src, size_t dest,
                                double lambda, size_t K, size_t* count_ptr) {
    // refer to: https://en.wikipedia.org/wiki/Yen%27s_algorithm
    auto g = graph;
    size_t count = 0;
    vector<GraphPath> A;
    A.reserve(K);
    auto [p1, p2] = biweight_dijkstra(g, src, dest, lambda);
    count++;
    A.emplace_back(move(p1));

    auto cmp = [&](const auto& lhs, const auto& rhs) {
        return graph.weighted_dist(lhs, lambda) >
               graph.weighted_dist(rhs, lambda);
    };
    priority_queue<GraphPath, deque<GraphPath>, decltype(cmp)> B(cmp);

    for (size_t k = 1; k < K; ++k) {
        for (size_t i = 0; i < A[k - 1].size() - 2; ++i) {
            auto spurNode = A[k - 1][i];
            GraphPath rootPath(A[k - 1].begin(), A[k - 1].begin() + i + 1);
            for (const auto& path : A) {
                bool flag = true;
                for (size_t j = 0; j <= i; ++j) {
                    if (path[j] != rootPath[j]) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    g.delete_edge(path[i], path[i + 1]);
                }
            }
            unordered_set<size_t> rootPathNode(rootPath.begin(),
                                               rootPath.end() - 1);
            g.delete_vertex(rootPathNode);

            auto [spurPath, _] = biweight_dijkstra(g, spurNode, dest, lambda);
            count++;
            if (!spurPath.empty()) {
                rootPath.insert(rootPath.end(), spurPath.begin() + 1,
                                spurPath.end());
                B.push(move(rootPath));
            }
            g = graph;
        }
        if (B.empty()) {
            break;
        }
        A.emplace_back(B.top());
        B.pop();
        while (!B.empty() && B.top() == A.back())
            B.pop();
    }

    if (count_ptr)
        *count_ptr += count;
    return A;
}
} // namespace BiLAD
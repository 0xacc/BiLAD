#define _USE_MATH_DEFINES
#include <bilad.h>
#include <cmath>

using namespace std;
namespace BiLAD {

// TODO: SCSP with p_1 & p_2

tuple<GraphPath, GraphPath, GraphPath, GraphPath, double, double, size_t,
      size_t>
new_bilad(TriWeightedGraph& graph, size_t src, size_t dest, weight_t &delta1,
          weight_t &delta2, SCSP scsp) {

    delta1 = delta1 * 2 + 1;
    delta2 = delta2 * 2 + 1;

    Graph g_c_d1 = graph.copy_c_d1();

    size_t dijkstra_count = 0;

    // S1
    auto [p_c, _] = biweight_dijkstra(g_c_d1, src, dest, 0);
    ++dijkstra_count;
    if (graph.d1_dist(p_c) < delta1 && graph.d2_dist(p_c) < delta2) {
        return {p_c, {}, {}, {}, 0, 0, 2, dijkstra_count};
    }

    // S2
    if (graph.d2_dist(p_c) < delta2) {
        graph.swap_d1_d2();
        swap(delta1, delta2);
    }
    g_c_d1 = graph.copy_c_d1();
    Graph g_c_d2 = graph.copy_c_d2(), g_d1_d2 = graph.copy_d1_d2();

    auto [p_2, _2] = biweight_dijkstra(g_c_d2, src, dest, 1, 0);
    ++dijkstra_count;
    if (graph.d2_dist(p_2) > delta2) {
        return {{}, {}, {}, {}, 0, 0, 0, dijkstra_count};
    }

    auto [p_plus, p_minus, lambda_2, count] = scsp(g_c_d2, src, dest, delta2);
    dijkstra_count += count;

    double alpha = (delta2 - graph.d2_dist(p_minus)) /
                   (double)(graph.d2_dist(p_plus) - graph.d2_dist(p_minus));
    double u_tilde =
        (1 - alpha) * graph.d1_dist(p_minus) + alpha * graph.d1_dist(p_plus);

    if (u_tilde <= delta1) {
        return {p_plus, p_minus, {}, {}, 0, lambda_2, 4, dijkstra_count};
    }

    double theta_1_low = 0, d1_plus = u_tilde,
           c_plus = (1 - alpha) * graph.c_dist(p_minus) +
                    alpha * graph.c_dist(p_plus);
    GraphPath p_minus_low = move(p_minus), p_plus_low = move(p_plus);

    // S3
    auto [p_1, _1] = biweight_dijkstra(g_c_d1, src, dest, 1, 0);
    ++dijkstra_count;
    if (graph.d1_dist(p_1) > delta1) {
        return {{}, {}, {}, {}, 0, 0, 0, dijkstra_count};
    }

    auto [p_plus2, p_minus2, lambda_22, count2] =
        scsp(g_d1_d2, src, dest, delta2);
    dijkstra_count += count2;

    double alpha2 = (delta2 - graph.d2_dist(p_minus2)) /
                    (double)(graph.d2_dist(p_plus2) - graph.d2_dist(p_minus2));
    double u_tilde2 = (1 - alpha2) * graph.d1_dist(p_minus2) +
                      alpha2 * graph.d1_dist(p_plus2);

    if (u_tilde2 >= delta1) {
        return {{}, {}, {}, {}, 0, 0, 0, dijkstra_count};
    }

    double theta_1_high = M_PI / 2, d1_minus = u_tilde2,
           c_minus = (1 - alpha) * graph.c_dist(p_minus2) +
                     alpha * graph.c_dist(p_plus2);
    GraphPath p_minus_high = move(p_minus2), p_plus_high = move(p_plus2);

    double lambda_1 = -1;

    const double gama = 0.05;

    // S4
    while (abs(theta_1_high - theta_1_low) > EPSILON) {
        lambda_1 = (c_plus - c_minus) / (d1_minus - d1_plus);
        double theta = atan(lambda_1);
        if (abs(theta - (theta_1_high + theta_1_low) / 2) >=
            (0.5 - gama) * (theta_1_high - theta_1_low)) {
            theta = (theta_1_high + theta_1_low) / 2;
            lambda_1 = tan(theta);
        }

        // S5
        auto [p_c_tilde, _3] = biweight_dijkstra(g_c_d1, src, dest, lambda_1);
        ++dijkstra_count;
        if (graph.d2_dist(p_c_tilde) < delta2) {
            if (graph.d1_dist(p_c_tilde) < delta1) {
                theta_1_high = theta;
                p_minus_high = p_c_tilde;
                p_plus_high = {};
            } else {
                theta_1_low = theta;
                p_minus_low = {};
                p_plus_low = p_c_tilde;
            }
            continue;
        }

        // S6
        Graph g_c_lambda_d1 = graph.copy_c_lambda_d1(lambda_1);
        auto [p_plus, p_minus, lambda, count] =
            scsp(g_c_lambda_d1, src, dest, delta2);
        dijkstra_count += count;
        double alpha = (delta2 - graph.d2_dist(p_minus)) /
                       (double)(graph.d2_dist(p_plus) - graph.d2_dist(p_minus));
        double u_tilde = (1 - alpha) * graph.d1_dist(p_minus) +
                         alpha * graph.d1_dist(p_plus);
        if (abs(u_tilde - delta1) < EPSILON) {
            return {p_plus,   p_minus, {}, {},
                    lambda_1, lambda,  1,  dijkstra_count};
        } else if (u_tilde < delta1) {
            theta_1_high = theta;
            d1_minus = u_tilde;
            c_minus = (1 - alpha) * graph.c_dist(p_minus) +
                      alpha * graph.c_dist(p_plus);
            p_minus_high = move(p_minus);
            p_plus_high = move(p_plus);
            lambda_2 = lambda;
        } else {
            theta_1_low = theta;
            d1_plus = u_tilde;
            c_plus = (1 - alpha) * graph.c_dist(p_minus) +
                     alpha * graph.c_dist(p_plus);
            p_minus_low = move(p_minus);
            p_plus_low = move(p_plus);
        }
    }
    return {p_plus_high, p_minus_high, p_plus_low, p_minus_low,
            lambda_1,    lambda_2,     3,          dijkstra_count};
}
} // namespace BiLAD
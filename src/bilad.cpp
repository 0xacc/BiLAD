#define _USE_MATH_DEFINES
#include <bilad.h>
#include <cmath>

using namespace std;
namespace BiLAD {
    std::tuple<GraphPath, GraphPath, double, bool, size_t>
    bilad(const Graph &graph, size_t src, size_t dest, weight_t delta) {
        size_t dijkstra_count = 0;
        const double decay = 0.05;
        bool flag = false;
        auto[_, p_plus]=biweight_dijkstra(graph, src, dest, 0);
        dijkstra_count++;
        weight_t delay_p_plus = graph.d_dist(p_plus);
        if (delay_p_plus <= delta) {
            return {p_plus, {}, -1, flag, dijkstra_count};
        }
        auto[p_minus, _1]=biweight_dijkstra(graph, src, dest, 1, 0);
        dijkstra_count++;
        weight_t delay_p_minus = graph.d_dist(p_minus);
        if (delay_p_minus > delta)return {{}, {}, -1, flag, dijkstra_count};
        if (delay_p_minus == delta)return {{}, p_minus, -1, flag, dijkstra_count};
        double theta_up = M_PI / 2., theta_down = 0;
        weight_t cost_p_plus, cost_p_minus;
        double lambda = -1;

        while (abs(theta_up - theta_down) > EPSILON) {
            cost_p_plus = graph.c_dist(p_plus);
            cost_p_minus = graph.c_dist(p_minus);
            delay_p_minus = graph.d_dist(p_minus);
            delay_p_plus = graph.d_dist(p_plus);

            lambda = static_cast<double>(cost_p_minus - cost_p_plus) / (delay_p_plus - delay_p_minus);
            double theta = atan(lambda);
            if (abs(theta - (theta_up + theta_down) / 2) > (0.5 - decay) * (theta_up - theta_down)) {
                theta = (theta_up + theta_down) / 2;
                lambda = tan(theta);
            }
            auto[p_c, p_d]=biweight_dijkstra(graph, src, dest, lambda);
            dijkstra_count++;
            weight_t delay_p_c = graph.d_dist(p_c), delay_p_d = graph.d_dist(p_d);
            if (delay_p_d <= delta && delta <= delay_p_c) {
                p_plus = p_c;
                p_minus = p_d;
                flag = delta - delay_p_minus == 0 || delay_p_plus - delta == 0;
                return {p_plus, p_minus, lambda, flag, dijkstra_count};
            } else if (delay_p_d > delta) {
                p_plus = p_d;
                theta_down = theta;
            } else if (delay_p_c < delta) {
                p_minus = p_c;
                theta_up = theta;
            }
        }
        return {p_plus, p_minus, lambda, true, dijkstra_count};
    }
}




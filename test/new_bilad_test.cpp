#include <bilad.h>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <random>
#include <string>

using namespace std;
using namespace BiLAD;

static TriWeightedGraph gen_graph(size_t n, double p = 0.99) {
    static auto e = default_random_engine();
    static uniform_real_distribution<double> link_type(0, 1);
    static uniform_int_distribution<uint32_t> long_delay(20, 30);
    static uniform_int_distribution<uint32_t> middle_delay(10, 20);
    static uniform_int_distribution<uint32_t> short_delay(1, 10);
    static uniform_int_distribution<uint32_t> cost(1, 30);

    TriWeightedGraph graph(n);

    GraphPath p1, p2;
    size_t count = 1;
    while (p1.empty() && p2.empty()) {
        // cout<<count++<<" tries..."<<endl;
        TriWeightedGraph result(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                if (link_type(e) < p)
                    continue;
                double r = link_type(e);
                weight_t link_cost = cost(e);
                weight_t link_delay;
                weight_t link_delay2;
                if (r < 0.75) {
                    link_delay = short_delay(e);
                } else if (r < 0.95) {
                    link_delay = middle_delay(e);
                } else {
                    link_delay = long_delay(e);
                }

                r = link_type(e);
                if (r < 0.75) {
                    link_delay2 = short_delay(e);
                } else if (r < 0.95) {
                    link_delay2 = middle_delay(e);
                } else {
                    link_delay2 = long_delay(e);
                }

                result.add_edge(i, j, link_cost, link_delay, link_delay2);
            }
        }
        auto path = biweight_dijkstra(result.copy_c_d1(), 0, n - 1, 0);
        p1 = move(get<0>(path));
        p2 = move(get<1>(path));
        graph = move(result);
    }

    return graph;
}

TEST(basic_test, all_feasible) {
    TriWeightedGraph m(3);
    m.add_edge(0, 1, 13, 0, 0);
    m.add_edge(0, 2, 10, 0, 0);
    m.add_edge(1, 2, 3, 0, 0);

    m.preprocess();

    weight_t delta1 = 1, delta2 = 1;

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 2, delta1, delta2, bind_bilad);

    EXPECT_EQ(p1, (GraphPath{0, 2}));
}

TEST(basic_test, no_d1_feasible) {
    TriWeightedGraph m(3);
    m.add_edge(0, 1, 13, 10, 0);
    m.add_edge(0, 2, 10, 10, 0);
    m.add_edge(1, 2, 3, 10, 0);

    m.preprocess();

    weight_t delta1 = 1, delta2 = 1;

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 2, delta1, delta2, bind_bilad);

    EXPECT_EQ(flag, 0);
}

TEST(basic_test, no_d2_feasible) {
    TriWeightedGraph m(3);
    m.add_edge(0, 1, 13, 0, 10);
    m.add_edge(0, 2, 10, 0, 10);
    m.add_edge(1, 2, 3, 0, 10);

    m.preprocess();

    weight_t delta1 = 1, delta2 = 1;

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 2, delta1, delta2, bind_bilad);

    EXPECT_EQ(flag, 0);
}

TEST(basic_test, feasible) {
    TriWeightedGraph m(3);
    m.add_edge(0, 1, 13, 10, 10);
    m.add_edge(0, 2, 10, 1, 1);
    m.add_edge(1, 2, 7, 1, 1);

    m.preprocess();

    weight_t delta1 = 5, delta2 = 5;

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 1, delta1, delta2, bind_bilad);

    EXPECT_EQ(flag, 4);
    EXPECT_EQ(p1, (GraphPath{0, 1}));
    EXPECT_EQ(p2, (GraphPath{0, 2, 1}));
}

TEST(basic_test, dual_optima) {
    TriWeightedGraph m(6);
    m.add_edge(0, 1, 2, 10, 1);
    m.add_edge(0, 3, 5, 1, 9);
    m.add_edge(0, 4, 20, 8, 8);
    m.add_edge(1, 2, 14, 5, 2);
    m.add_edge(1, 5, 3, 15, 24);
    m.add_edge(2, 5, 14, 5, 2);
    m.add_edge(3, 4, 5, 1, 9);
    m.add_edge(4, 5, 20, 2, 2);

    m.preprocess();

    weight_t delta1 = 15, delta2 = 15;

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 5, delta1, delta2, bind_bilad);

    EXPECT_EQ(flag, 3);
    EXPECT_EQ(p1, (GraphPath{0, 3, 4, 5}));
    EXPECT_EQ(p2, (GraphPath{0, 1, 2, 5}));
    EXPECT_EQ(p3, (GraphPath{0, 1, 5}));
    EXPECT_EQ(p4, (GraphPath{0, 1, 2, 5}));
}

TEST(basic_test, dual_optima_infeasible) {
    TriWeightedGraph m(6);
    m.add_edge(0, 1, 2, 10, 1);
    m.add_edge(0, 3, 5, 1, 9);
    // m.add_edge(0, 4, 20, 8, 8);
    m.add_edge(1, 2, 14, 5, 2);
    m.add_edge(1, 5, 3, 15, 24);
    m.add_edge(2, 5, 14, 5, 2);
    m.add_edge(3, 4, 5, 1, 9);
    m.add_edge(4, 5, 20, 2, 2);

    m.preprocess();

    weight_t delta1 = 15, delta2 = 15;

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 5, delta1, delta2, bind_bilad);

    EXPECT_EQ(flag, 3);
    EXPECT_EQ(p1, (GraphPath{0, 3, 4, 5}));
    EXPECT_EQ(p2, (GraphPath{0, 1, 2, 5}));
    EXPECT_EQ(p3, (GraphPath{0, 1, 5}));
    EXPECT_EQ(p4, (GraphPath{0, 1, 2, 5}));
}

TEST(radom_test, small_graph) {
    for (int i = 0; i < 100; ++i) {
        TriWeightedGraph m = gen_graph(500);

        auto [p_d2_1, p_d2_2] =
            biweight_dijkstra(m.copy_d1_d2(), 0, m.size() - 1, 1, 0);
        EXPECT_NE(p_d2_1, (GraphPath{}));
        auto delta2 = 2 * m.d2_dist(p_d2_1);
        auto [p_d1, p_d2, lambda, _] =
            bind_bilad(m.copy_d1_d2(), 0, m.size() - 1, delta2);
        EXPECT_NE(p_d1, (GraphPath{}));
        auto delta1 = 2 * m.d1_dist(p_d1);

        m.preprocess();

        auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
            new_bilad(m, 0, m.size() - 1, delta1, delta2, bind_bilad);

        cout << "id: " << i << endl;
        cout << "flag: " << flag << endl;
        cout << "count: " << count << endl;

        if (flag == 3) {
            if (!p1.empty() && !p2.empty()) {
                EXPECT_GT(m.d2_dist(p1), delta2);
                EXPECT_LT(m.d2_dist(p2), delta2);
            }
            if (!p3.empty() && !p4.empty()) {
                EXPECT_GT(m.d2_dist(p3), delta2);
                EXPECT_LT(m.d2_dist(p4), delta2);
            }
        }
    }
}

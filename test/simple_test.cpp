#include <gtest/gtest.h>
#include<bilad.h>
#include<fstream>
#include<string>
#include<iostream>

using namespace BiLAD;

TEST(simple_test, sigle_weight_direct_path) {
    Graph m(3);
    m.add_edge(0, 1, 13, 0);
    m.add_edge(0, 2, 10, 0);
    m.add_edge(1, 2, 3, 0);
    auto[p1, p2]=biweight_dijkstra(m, 0, 2, 0);
    EXPECT_EQ(p1, (GraphPath{0, 2}));
    EXPECT_EQ(p2, (GraphPath{0, 2}));
    EXPECT_EQ(m.d_dist(p1), 0);
    EXPECT_EQ(m.d_dist(p2), 0);
    EXPECT_EQ(m.c_dist(p1), 10);
    EXPECT_EQ(m.c_dist(p2), 10);
    EXPECT_EQ(m.weighted_dist(p1, 1), 10);
    EXPECT_EQ(m.weighted_dist(p2, 1), 10);
}

TEST(simple_test, sigle_weight_intermediate_path) {
    Graph m(3);
    m.add_edge(0, 1, 3, 0);
    m.add_edge(0, 2, 10, 0);
    m.add_edge(1, 2, 3, 0);
    auto[p1, p2]=biweight_dijkstra(m, 0, 2, 0);
    EXPECT_EQ(p1, (GraphPath{0, 1, 2}));
    EXPECT_EQ(p2, (GraphPath{0, 1, 2}));
    EXPECT_EQ(m.d_dist(p1), 0);
    EXPECT_EQ(m.d_dist(p2), 0);
    EXPECT_EQ(m.c_dist(p1), 6);
    EXPECT_EQ(m.c_dist(p2), 6);
    EXPECT_EQ(m.weighted_dist(p1, 1), 6);
    EXPECT_EQ(m.weighted_dist(p2, 1), 6);
}

TEST(simple_test, double_weight_test_1) {
    Graph m(5);
    m.add_edge(0, 4, 10, 10);
    m.add_edge(0, 2, 1, 1);
    m.add_edge(2, 4, 1, 1);
    m.add_edge(0, 1, 1, 2);
    m.add_edge(1, 4, 2, 1);
    m.add_edge(0, 3, 3, 3);
    m.add_edge(3, 4, 1, 1);
    m.add_edge(1, 2, 1, 1);
    m.add_edge(2, 3, 1, 1);
    auto[p1, p2]=biweight_dijkstra(m, 0, 4, 1);
    EXPECT_EQ(p1, (GraphPath{0, 2, 4}));
    EXPECT_EQ(p2, (GraphPath{0, 2, 4}));
    EXPECT_EQ(m.d_dist(p1), 2);
    EXPECT_EQ(m.d_dist(p2), 2);
    EXPECT_EQ(m.c_dist(p1), 2);
    EXPECT_EQ(m.c_dist(p2), 2);
    EXPECT_EQ(m.weighted_dist(p1, 1), 4);
    EXPECT_EQ(m.weighted_dist(p2, 1), 4);
}

TEST(simple_test, double_weight_test_2) {
    Graph m(5);
    m.add_edge(0, 4, 10, 10);
    m.add_edge(0, 2, 1, 1);
    m.add_edge(2, 4, 1, 1);
    m.add_edge(0, 1, 2, 0);
    m.add_edge(1, 4, 2, 0);
    m.add_edge(0, 3, 0, 2);
    m.add_edge(3, 4, 0, 2);
    m.add_edge(1, 2, 1, 1);
    m.add_edge(2, 3, 1, 1);
    auto[p1, p2]=biweight_dijkstra(m, 0, 4, 1);
    EXPECT_EQ(p1, (GraphPath{0, 3, 4}));
    EXPECT_EQ(p2, (GraphPath{0, 1, 4}));
    EXPECT_EQ(m.d_dist(p1), 4);
    EXPECT_EQ(m.d_dist(p2), 0);
    EXPECT_EQ(m.c_dist(p1), 0);
    EXPECT_EQ(m.c_dist(p2), 4);
    EXPECT_EQ(m.weighted_dist(p1, 1), 4);
    EXPECT_EQ(m.weighted_dist(p2, 1), 4);
}

TEST(simple_test, double_weight_test_3) {
    {
        Graph m(6);
        m.add_edge(0, 2, 2, 2);
        m.add_edge(2, 3, 2, 2);
        m.add_edge(3, 5, 2, 2);
        m.add_edge(0, 3, 5, 5);
        m.add_edge(2, 5, 5, 5);

        auto[p1, p2]=biweight_dijkstra(m, 0, 5, 1);
        EXPECT_EQ(p1, (GraphPath{0, 2, 3, 5}));
        EXPECT_EQ(p2, (GraphPath{0, 2, 3, 5}));
        EXPECT_EQ(m.d_dist(p1), 6);
        EXPECT_EQ(m.d_dist(p2), 6);
        EXPECT_EQ(m.c_dist(p1), 6);
        EXPECT_EQ(m.c_dist(p2), 6);
        EXPECT_EQ(m.weighted_dist(p1, 1), 12);
        EXPECT_EQ(m.weighted_dist(p2, 1), 12);
    }
    {
        Graph m(6);
        m.add_edge(0, 2, 1, 1);
        m.add_edge(2, 3, 1, 1);
        m.add_edge(3, 5, 1, 1);
        m.add_edge(2, 5, 3, 1);

        auto[p1, p2]=biweight_dijkstra(m, 0, 5, 1);
        EXPECT_EQ(p1, (GraphPath{0, 2, 3, 5}));
        EXPECT_EQ(p2, (GraphPath{0, 2, 5}));
        EXPECT_EQ(m.d_dist(p1), 3);
        EXPECT_EQ(m.d_dist(p2), 2);
        EXPECT_EQ(m.c_dist(p1), 3);
        EXPECT_EQ(m.c_dist(p2), 4);
        EXPECT_EQ(m.weighted_dist(p1, 1), 6);
        EXPECT_EQ(m.weighted_dist(p2, 1), 6);
    }
}

TEST(simple_test, double_weight_external) {
    for (int i = 1; i <= 10; ++i) {
        std::string filename = "data/Id_" + std::to_string(i) + ".txt";
        std::fstream f(filename);
        if (!f) {
            std::cerr << "failed to open " << filename << std::endl;
            EXPECT_EQ(static_cast<bool>(f), true);
        }
        std::string b;
        f >> b >> b >> b >> b >> b;
        Graph m(20);
        int id, src, dest, c, d;
        while (f >> id >> src >> dest >> c >> d) {
            m.add_edge(src, dest, c, d);
        }

        {
            auto[p1, p2]=biweight_dijkstra(m, 3, 19, 0);
            for (auto v:p1) {
                std::cout << v << "->";
            }
            std::cout << '[' << "cost: " << m.c_dist(p1) << ", delay: " << m.d_dist(p1) << ']';
            std::cout << '\n';
            for (auto v:p2) {
                std::cout << v << "->";
            }
            std::cout << '[' << "cost: " << m.c_dist(p2) << ", delay: " << m.d_dist(p2) << ']';
            std::cout << '\n';
        }
        {
            auto[p1, p2]=biweight_dijkstra(m, 3, 19, 1, 0);
            for (auto v:p1) {
                std::cout << v << "->";
            }
            std::cout << '[' << "cost: " << m.c_dist(p1) << ", delay: " << m.d_dist(p1) << ']';
            std::cout << '\n';
            for (auto v:p2) {
                std::cout << v << "->";
            }
            std::cout << '[' << "cost: " << m.c_dist(p2) << ", delay: " << m.d_dist(p2) << ']';
            std::cout << '\n';
        }
        auto[p_plus, p_minus, lambda, flag, count]=bilad(m, 3, 19, 24);
        std::cout << count << " times dijkstra calls" << '\n';
        if (!flag) {
            if (p_plus.empty() && p_minus.empty()) {
                std::cout << "infeasible!";
            } else {
                if (!p_plus.empty() && !p_minus.empty()) {
                    for (auto v:p_plus)std::cout << v << "->";
                    std::cout << '[' << "cost: " << m.c_dist(p_plus) << ", delay: " << m.d_dist(p_plus) << ']';
                    std::cout << '\n';
                    for (auto v:p_minus)std::cout << v << "->";
                    std::cout << '[' << "cost: " << m.c_dist(p_minus) << ", delay: " << m.d_dist(p_minus) << ']';
                } else if (!p_plus.empty()) {
                    for (auto v:p_plus)std::cout << v << "->";
                    std::cout << '[' << "cost: " << m.c_dist(p_plus) << ", delay: " << m.d_dist(p_plus) << ']';
                } else if (!p_minus.empty()) {
                    for (auto v:p_minus)std::cout << v << "->";
                    std::cout << '[' << "cost: " << m.c_dist(p_minus) << ", delay: " << m.d_dist(p_minus) << ']';
                }
            }
        } else {
            std::cout << "gap!\n";
            for (auto v:p_plus)std::cout << v << "->";
            std::cout << '[' << "cost: " << m.c_dist(p_plus) << ", delay: " << m.d_dist(p_plus) << ']';
            std::cout << '\n';
            for (auto v:p_minus)std::cout << v << "->";
            std::cout << '[' << "cost: " << m.c_dist(p_minus) << ", delay: " << m.d_dist(p_minus) << ']';
        }
        std::cout << std::endl << std::endl;
    }
}

TEST(simple_test, exact_bilad) {
    for (int i = 1; i <= 10; ++i) {
        std::string filename = "data/Id_" + std::to_string(i) + ".txt";
        std::fstream f(filename);
        if (!f) {
            std::cerr << "failed to open " << filename << std::endl;
            EXPECT_EQ(static_cast<bool>(f), true);
        }
        std::string b;
        f >> b >> b >> b >> b >> b;
        Graph m(20);
        int id, src, dest, c, d;
        while (f >> id >> src >> dest >> c >> d) {
            m.add_directed_edge(src, dest, c, d);
        }

        {
            auto[p1, p2]=biweight_dijkstra(m, 3, 19, 0);
            for (auto v:p1) {
                std::cout << v << "->";
            }
            std::cout << '[' << "cost: " << m.c_dist(p1) << ", delay: " << m.d_dist(p1) << ']';
            std::cout << '\n';
            for (auto v:p2) {
                std::cout << v << "->";
            }
            std::cout << '[' << "cost: " << m.c_dist(p2) << ", delay: " << m.d_dist(p2) << ']';
            std::cout << '\n';
        }
        {
            auto[p1, p2]=biweight_dijkstra(m, 3, 19, 1, 0);
            for (auto v:p1) {
                std::cout << v << "->";
            }
            std::cout << '[' << "cost: " << m.c_dist(p1) << ", delay: " << m.d_dist(p1) << ']';
            std::cout << '\n';
            for (auto v:p2) {
                std::cout << v << "->";
            }
            std::cout << '[' << "cost: " << m.c_dist(p2) << ", delay: " << m.d_dist(p2) << ']';
            std::cout << '\n';
        }
        auto[p, count]=exact_bilad(m, 3, 19, 24);
        std::cout << count << " times dijkstra calls" << '\n';
        for (auto v:p)std::cout << v << "->";
        std::cout << '[' << "cost: " << m.c_dist(p) << ", delay: " << m.d_dist(p) << ']';
        std::cout << '\n';
        std::cout << std::endl;
    }
}

TEST(simple_test, yen_algorithm) {
    Graph m(6);
    m.add_edge(0, 1, 3, 0);
    m.add_edge(1, 2, 4, 0);
    m.add_edge(0, 3, 2, 0);
    m.add_edge(1, 3, 1, 0);
    m.add_edge(2, 3, 2, 0);
    m.add_edge(2, 4, 2, 0);
    m.add_edge(3, 4, 3, 0);
    m.add_edge(2, 5, 1, 0);
    m.add_edge(4, 5, 2, 0);

    auto v = yen_algorithm(m, 0, 5, 0, 4);

    for (const auto &path:v) {
        for (const auto &node:path) {
            std::cout << node << "->";
        }
        std::cout << std::endl;
    }
    EXPECT_EQ(v[0], (GraphPath{0, 3, 2, 5}));
    EXPECT_EQ(m.weighted_dist(v[1], 1), 7);
    EXPECT_EQ(m.weighted_dist(v[2], 1), 7);
    EXPECT_EQ(v[3], (GraphPath{0, 1, 2, 5}));
}

TEST(simple_test, yen_algorithm_directed) {
    Graph m(6);
    m.add_directed_edge(0, 2, 3, 0);
    m.add_directed_edge(2, 4, 4, 0);
    m.add_directed_edge(0, 1, 2, 0);
    m.add_directed_edge(1, 2, 1, 0);
    m.add_directed_edge(1, 4, 2, 0);
    m.add_directed_edge(1, 3, 3, 0);
    m.add_directed_edge(4, 3, 2, 0);
    m.add_directed_edge(3, 5, 2, 0);
    m.add_directed_edge(4, 5, 1, 0);

    auto v = yen_algorithm(m, 0, 5, 0, 4);

    for (const auto &path:v) {
        for (const auto &node:path) {
            std::cout << node << "->";
        }
        std::cout << std::endl;
    }
    EXPECT_EQ(v[0], (GraphPath{0, 1, 4, 5}));
    EXPECT_EQ(v[1], (GraphPath{0, 1, 3, 5}));
    EXPECT_EQ(v[2], (GraphPath{0, 2, 4, 5}));
    EXPECT_GE(m.c_dist(v[3]), 8);
}
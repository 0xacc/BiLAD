#include <bilad.h>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <string>

using namespace BiLAD;

TEST(basic_test, all_feasible) {
    TriWeightedGraph m(3);
    m.add_edge(0, 1, 13, 0, 0);
    m.add_edge(0, 2, 10, 0, 0);
    m.add_edge(1, 2, 3, 0, 0);

    m.preprocess();

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 2, 1, 1, bind_bilad);

    EXPECT_EQ(p1, (GraphPath{0, 2}));
}

TEST(basic_test, no_d1_feasible) {
    TriWeightedGraph m(3);
    m.add_edge(0, 1, 13, 10, 0);
    m.add_edge(0, 2, 10, 10, 0);
    m.add_edge(1, 2, 3, 10, 0);

    m.preprocess();

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 2, 1, 1, bind_bilad);

    EXPECT_EQ(flag, 0);
}

TEST(basic_test, no_d2_feasible) {
    TriWeightedGraph m(3);
    m.add_edge(0, 1, 13, 0, 10);
    m.add_edge(0, 2, 10, 0, 10);
    m.add_edge(1, 2, 3, 0, 10);

    m.preprocess();

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 2, 1, 1, bind_bilad);

    EXPECT_EQ(flag, 0);
}

TEST(basic_test, feasible) {
    TriWeightedGraph m(3);
    m.add_edge(0, 1, 13, 10, 10);
    m.add_edge(0, 2, 10, 1, 1);
    m.add_edge(1, 2, 7, 1, 1);

    m.preprocess();

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 1, 5, 5, bind_bilad);

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

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 5, 15, 15, bind_bilad);

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

    auto [p1, p2, p3, p4, lambda1, lambda2, flag, count] =
        new_bilad(m, 0, 5, 15, 15, bind_bilad);

    EXPECT_EQ(flag, 3);
    EXPECT_EQ(p1, (GraphPath{0, 3, 4, 5}));
    EXPECT_EQ(p2, (GraphPath{0, 1, 2, 5}));
    EXPECT_EQ(p3, (GraphPath{0, 1, 5}));
    EXPECT_EQ(p4, (GraphPath{0, 1, 2, 5}));
}
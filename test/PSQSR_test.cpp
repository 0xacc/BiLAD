#include <gtest/gtest.h>
#include<bilad.h>
#include<fstream>
#include<string>
#include<iostream>

using namespace BiLAD;

// TODO: MAKE THE GAP

TEST(m_larac, sigle_weight_direct_path) {
    Graph m(3);
    m.add_edge(0, 1, 13, 0);
    m.add_edge(0, 2, 10, 0);
    m.add_edge(1, 2, 3, 0);
    auto [p1, p2,count]=PSQSR(m,0,2,0);
    EXPECT_EQ(p1, (GraphPath{0, 2}));
    EXPECT_EQ(p2, (GraphPath{}));
    EXPECT_EQ(m.d_dist(p1), 0);
    EXPECT_EQ(m.c_dist(p1), 10);
    EXPECT_EQ(m.weighted_dist(p1, 1), 10);
}

TEST(m_larac, sigle_weight_intermediate_path) {
    Graph m(3);
    m.add_edge(0, 1, 3, 0);
    m.add_edge(0, 2, 10, 0);
    m.add_edge(1, 2, 3, 0);
    auto[p1, p2, count]=PSQSR(m, 0, 2, 0);
    EXPECT_EQ(p1, (GraphPath{0, 1, 2}));
    EXPECT_EQ(p2, (GraphPath{}));
    EXPECT_EQ(m.d_dist(p1), 0);
    EXPECT_EQ(m.c_dist(p1), 6);
    EXPECT_EQ(m.weighted_dist(p1, 1), 6);
}

TEST(m_larac, double_weight_test_1) {
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
    auto[p1, p2, count]=PSQSR(m, 0, 4, 1);
    EXPECT_EQ(p1, (GraphPath{0, 2, 4}));
    EXPECT_EQ(m.d_dist(p1), 2);
    EXPECT_EQ(m.c_dist(p1), 2);
    EXPECT_EQ(m.weighted_dist(p1, 1), 4);
}
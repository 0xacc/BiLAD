#include <gtest/gtest.h>
#include <bilad.h>
#include <random>
#include <iostream>
#include <chrono>

using namespace std;
using namespace BiLAD;

Graph gen_graph(size_t n){
    static auto e=default_random_engine();
    static uniform_real_distribution<double> link_type(0,1);
    static uniform_int_distribution<uint32_t> long_delay(20,30);
    static uniform_int_distribution<uint32_t> middle_delay(10,20);
    static uniform_int_distribution<uint32_t> short_delay(1,10);
    static uniform_int_distribution<uint32_t> cost(1,30);


    Graph result(n);

    for(size_t i=0;i<n;++i){
        for(size_t j=i+1;j<n;++j){
            //if(link_type(e)<0.5)continue;
            double r=link_type(e);
            weight_t link_cost=cost(e);
            weight_t link_delay;
            if(r<0.75){
                link_delay=short_delay(e);
            }else if(r<0.95){
                link_delay=middle_delay(e);
            }else{
                link_delay=long_delay(e);
            }

            result.add_edge(i,j,link_cost,link_delay);
        }
    }

    return result;
}

TEST(random_test,small_net){
    Graph g=gen_graph(200);
    auto [p1,p2]=biweight_dijkstra(g,0,9,1,0);
    ASSERT_EQ(g.d_dist(p1),g.d_dist(p2));
    weight_t delta=g.d_dist(p1);
    delta*=2;

    auto t1=chrono::high_resolution_clock::now();
    auto [p_plus,p_minus,lambda,flag,count]=bilad(g, 0, 9, delta);
    auto t2=chrono::high_resolution_clock::now();

    auto duration=chrono::duration_cast<chrono::milliseconds>(t2-t1);

    cout<<"time: "<<duration.count()<<" ms"<<'\n';
    if(flag){
        cout<<"gap!"<<endl;
    }
    cout<<"p+ length: "<<p_plus.size()<<'\n';
    cout<<"p+ cost: "<<g.c_dist(p_plus)<<'\n';
    cout<<"p+ delay: "<<g.d_dist(p_plus)<<'\n';
    cout<<"p- length: "<<p_minus.size()<<'\n';
    cout<<"p- cost: "<<g.c_dist(p_minus)<<'\n';
    cout<<"p- delay: "<<g.d_dist(p_minus)<<'\n';
    cout<<"dijkstra calls: "<<count<<endl;
}

TEST(random_test,small_net_exact){
    Graph g=gen_graph(200);
    auto [p1,p2]=biweight_dijkstra(g,0,9,1,0);
    ASSERT_EQ(g.d_dist(p1),g.d_dist(p2));
    weight_t delta=g.d_dist(p1);
    delta*=2;

    auto t1=chrono::high_resolution_clock::now();
    auto [p,count]=exact_bilad(g, 0, 9, delta);
    auto t2=chrono::high_resolution_clock::now();

    auto duration=chrono::duration_cast<chrono::milliseconds>(t2-t1);

    cout<<"time: "<<duration.count()<<" ms"<<'\n';
    cout<<"p length: "<<p.size()<<'\n';
    cout<<"p cost: "<<g.c_dist(p)<<'\n';
    cout<<"p delay: "<<g.d_dist(p)<<'\n';
    cout<<"dijkstra calls: "<<count<<endl;
}
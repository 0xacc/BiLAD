#include <gtest/gtest.h>
#include <bilad.h>
#include <random>
#include <iostream>
#include <chrono>
#include <fstream>

using namespace std;
using namespace BiLAD;

Graph gen_graph(size_t n, double p=0.99){
    static auto e=default_random_engine();
    static uniform_real_distribution<double> link_type(0,1);
    static uniform_int_distribution<uint32_t> long_delay(20,30);
    static uniform_int_distribution<uint32_t> middle_delay(10,20);
    static uniform_int_distribution<uint32_t> short_delay(1,10);
    static uniform_int_distribution<uint32_t> cost(1,30);


    Graph graph(n);

    GraphPath p1,p2;
    size_t count=1;
    while(p1.empty()&&p2.empty()){
        //cout<<count++<<" tries..."<<endl;
        Graph result(n);
        for(size_t i=0;i<n;++i){
            for(size_t j=i+1;j<n;++j){
                if(link_type(e)<p)continue;
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
        auto path=biweight_dijkstra(result,0,n-1,0);
        p1=move(get<0>(path));
        p2=move(get<1>(path));
        graph=move(result);
    }

    return graph;
}

TEST(random_test,small_net){
    size_t n=2000;
    Graph g=gen_graph(n);// 2000
    auto [p1,p2]=biweight_dijkstra(g,0,n-1,1,0);
    ASSERT_EQ(g.d_dist(p1),g.d_dist(p2));
    weight_t delta=g.d_dist(p1);
    delta*=2;

    auto t1=chrono::high_resolution_clock::now();
    auto [p_plus,p_minus,lambda,flag,count]=bilad(g, 0, n-1, delta);
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
    size_t n=2000;
    Graph g=gen_graph(n,0.985);
    auto [p1,p2]=biweight_dijkstra(g,0,n-1,1,0);
    ASSERT_EQ(g.d_dist(p1),g.d_dist(p2));
    weight_t delta=g.d_dist(p1);
    delta*=2;

    auto t1=chrono::high_resolution_clock::now();
    auto [p,count]=exact_bilad(g, 0, n-1, delta);
    auto t2=chrono::high_resolution_clock::now();

    auto duration=chrono::duration_cast<chrono::milliseconds>(t2-t1);

    cout<<"time: "<<duration.count()<<" ms"<<'\n';
    cout<<"p length: "<<p.size()<<'\n';
    cout<<"p cost: "<<g.c_dist(p)<<'\n';
    cout<<"p delay: "<<g.d_dist(p)<<'\n';
    cout<<"dijkstra calls: "<<count<<endl;
}

TEST(random_test,plotting){
    fstream bilad_file("bilad.csv",fstream::out),exact_bilad_file("exact_bilad.csv",fstream::out);
    bilad_file<<"nodes"<<','<<"dijkstra"<<','<<"time"<<'\n';
    exact_bilad_file<<"nodes"<<','<<"dijkstra"<<','<<"time"<<'\n';
    for(size_t nodes=20;nodes<=2000;nodes+=20){
        double bilad_time=0;
        double exact_bilad_time=0;
        double bilad_dijkstra=0;
        double exact_bilad_dijkstra=0;

        constexpr size_t reapt_times=150;
        for(size_t i=0;i<reapt_times;++i){
            constexpr double degree=30;
            Graph g=gen_graph(nodes,1-degree/nodes);
            auto [p1,p2]=biweight_dijkstra(g,0,nodes-1,1,0);
            ASSERT_EQ(g.d_dist(p1),g.d_dist(p2));
            weight_t delta=g.d_dist(p1);
            delta*=2;

            auto t1=chrono::high_resolution_clock::now();
            auto [p_plus,p_minus,lambda,flag,bilad_count]=bilad(g, 0, nodes-1, delta);
            auto t2=chrono::high_resolution_clock::now();
            bilad_time+=chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
            bilad_dijkstra+=bilad_count;

            auto t3=chrono::high_resolution_clock::now();
            auto [p,exact_count]=exact_bilad(g, 0, nodes-1, delta);
            auto t4=chrono::high_resolution_clock::now();
            exact_bilad_time+=chrono::duration_cast<chrono::milliseconds>(t4-t3).count();
            exact_bilad_dijkstra+=exact_count;
        }
        bilad_file<<nodes<<','<<bilad_dijkstra/reapt_times<<','<<bilad_time/reapt_times<<'\n';
        exact_bilad_file<<nodes<<','<<exact_bilad_dijkstra/reapt_times<<','<<exact_bilad_time/reapt_times<<'\n';
    }
}

TEST(random_test,plotting_degree){
    fstream bilad_file("bilad_degree.csv",fstream::out),exact_bilad_file("exact_bilad_degree.csv",fstream::out);
    bilad_file<<"degree"<<','<<"dijkstra"<<','<<"time"<<'\n';
    exact_bilad_file<<"degree"<<','<<"dijkstra"<<','<<"time"<<'\n';

    constexpr size_t nodes=100;
    for(size_t degree=1;degree<=100;degree+=1){
        double bilad_time=0;
        double exact_bilad_time=0;
        double bilad_dijkstra=0;
        double exact_bilad_dijkstra=0;

        constexpr size_t reapt_times=150;
        for(size_t i=0;i<reapt_times;++i){
            Graph g=gen_graph(nodes,1-degree/static_cast<double>(nodes));
            auto [p1,p2]=biweight_dijkstra(g,0,nodes-1,1,0);
            ASSERT_EQ(g.d_dist(p1),g.d_dist(p2));
            weight_t delta=g.d_dist(p1);
            delta*=2;

            auto t1=chrono::high_resolution_clock::now();
            auto [p_plus,p_minus,lambda,flag,bilad_count]=bilad(g, 0, nodes-1, delta);
            auto t2=chrono::high_resolution_clock::now();
            bilad_time+=chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
            bilad_dijkstra+=bilad_count;

            auto t3=chrono::high_resolution_clock::now();
            auto [p,exact_count]=exact_bilad(g, 0, nodes-1, delta);
            auto t4=chrono::high_resolution_clock::now();
            exact_bilad_time+=chrono::duration_cast<chrono::milliseconds>(t4-t3).count();
            exact_bilad_dijkstra+=exact_count;
        }
        bilad_file<<degree<<','<<bilad_dijkstra/reapt_times<<','<<bilad_time/reapt_times<<'\n';
        exact_bilad_file<<degree<<','<<exact_bilad_dijkstra/reapt_times<<','<<exact_bilad_time/reapt_times<<'\n';
    }
}
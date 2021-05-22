#ifndef BILAD_BILAD_H
#define BILAD_BILAD_H

#include <vector>
#include <tuple>
#include <utility>
#include <cstdint>
#include <algorithm>

namespace BiLAD {
    const double EPSILON=1e-5;
    using weight_t = unsigned int;
    using GraphPath = std::vector<size_t>;

    class Graph{
    private:
        using AdjTable = std::vector<std::vector<std::pair<size_t, std::pair<weight_t, weight_t>>>>;
        AdjTable adj;
    public:
        explicit Graph(size_t n):adj(n){}
        void add_edge(size_t i, size_t j, weight_t c, weight_t d){
            adj[i].emplace_back(std::make_pair(j, std::make_pair(c, d)));
            adj[j].emplace_back(std::make_pair(i, std::make_pair(c, d)));
        }
        auto& operator[](size_t i)const{
            return adj[i];
        }
        size_t size()const{
            return adj.size();
        }
        weight_t d_dist(const GraphPath &path)const{
            weight_t dist=0;
            for(size_t i=1;i<path.size();++i){
                const auto &edges=adj[path[i-1]];
                const auto &target=*std::find_if(edges.begin(), edges.end(),[&](auto v){return v.first==path[i];});
                dist+=target.second.second;
            }
            return dist;
        }
        weight_t c_dist(const GraphPath &path)const{
            weight_t dist=0;
            for(size_t i=1;i<path.size();++i){
                const auto &edges=adj[path[i-1]];
                const auto &target=*std::find_if(edges.begin(), edges.end(),[&](auto v){return v.first==path[i];});
                dist+=target.second.first;
            }
            return dist;
        }
    };

    std::tuple<GraphPath, GraphPath> biweight_dijkstra(const Graph &graph, size_t src, size_t dest, double lambda, double gama=1);

    std::tuple<GraphPath, GraphPath, double, bool,size_t> bilad(const Graph &graph, size_t src, size_t dest, weight_t delta);

    GraphPath exact_bilad(const Graph &graph, size_t src, size_t dest, double delta);

}
#endif //BILAD_BILAD_H

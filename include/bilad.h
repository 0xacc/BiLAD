#ifndef BILAD_BILAD_H
#define BILAD_BILAD_H

#include <vector>
#include <tuple>
#include <utility>
#include <cstdint>

namespace BiLAD {
    using weight_t = unsigned int;
    using GraphPath = std::vector<size_t>;

    class Graph{
    private:
        using AdjTable = std::vector<std::vector<std::pair<size_t, std::pair<weight_t, weight_t>>>>;
        AdjTable adj;
    public:
        Graph(size_t n):adj(n){}
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
    };

    std::tuple<GraphPath, GraphPath> biweight_dijkstra(const Graph &graph, size_t src, size_t dest, double lambda, double gama=1);

    std::tuple<GraphPath, bool> bilad(const Graph &graph);

    GraphPath exact_bilad(const Graph &graph);

}
#endif //BILAD_BILAD_H

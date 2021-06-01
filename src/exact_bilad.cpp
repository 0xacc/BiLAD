#include <bilad.h>
#include <vector>
#include <unordered_set>
#include <queue>

using namespace std;

namespace BiLAD{
GraphPath exact_bilad(const Graph &graph, size_t src, size_t dest, double delta) {
    return GraphPath();
}

std::vector<GraphPath> yen_algorithm(const Graph &graph, size_t src, size_t dest, double lambda, size_t K) {
    // refer to: https://en.wikipedia.org/wiki/Yen%27s_algorithm
    auto g=graph;
    vector<GraphPath> A;
    A.reserve(K);
    auto [p1,p2]=biweight_dijkstra(g,src,dest,lambda);
    A.emplace_back(move(p1));

    auto cmp=[&](const auto &lhs,const auto &rhs){
        return graph.weighted_dist(lhs,lambda)>graph.weighted_dist(rhs,lambda);
    };
    priority_queue<GraphPath,deque<GraphPath>, decltype(cmp)> B(cmp);

    for(size_t k=1;k<K;++k){
        for(size_t i=0;i<A[k-1].size()-2;++i){
            auto spurNode=A[k-1][i];
            GraphPath rootPath(A[k-1].begin(),A[k-1].begin()+i+1);
            for(const auto &path:A){
                bool flag= true;
                for(size_t j=0;j<=i;++j){
                    if(path[i]!=rootPath[i]){
                        flag=false;
                        break;
                    }
                }
                if(flag){
                    g.delete_edge(path[i],path[i+1]);
                }
            }
            unordered_set<size_t> rootPathNode(rootPath.begin(), rootPath.end()-1);
            g.delete_vertex(rootPathNode);

            auto [spurPath,_]=biweight_dijkstra(g,spurNode,dest,lambda);
            if(!spurPath.empty()) {
                rootPath.insert(rootPath.end(), spurPath.begin() + 1, spurPath.end());
                B.push(move(rootPath));
            }

            g=graph;
        }
        if(B.empty()){
            break;
        }
        A.emplace_back(B.top());
        B.pop();
    }

    return A;
    }
}
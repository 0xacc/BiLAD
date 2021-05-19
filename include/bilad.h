#ifndef BILAD_BILAD_H
#define BILAD_BILAD_H
#include <vector>
#include <tuple>

using AdjMat=std::vector<std::vector<double>>;
using GraphPath=std::vector<size_t>;
std::tuple<GraphPath,GraphPath> biweight_dijkstra(const AdjMat &m1, const AdjMat &m2,size_t src,size_t dst,double lambda);
std::tuple<GraphPath,bool> bilad(const AdjMat &m1,const AdjMat &m2);
GraphPath exact_bilad(const AdjMat &m1,const AdjMat &m2);

#endif //BILAD_BILAD_H

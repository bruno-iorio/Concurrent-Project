#ifndef GRAPH_HPP
#define GRAPH_HPP
#pragma once

#include <vector>
#include <utility>   
#include <string>

class Graph {
public:
    int n = 0;                                           
    int maxDist = 0;                                            
    std::vector<std::vector<std::pair<int,double>>> adj_lists;     

    void parse_graph(const std::string& filename);            
};

struct customCompare {
    bool operator()(const std::pair<int,double>& a,
                    const std::pair<int,double>& b) const noexcept {
        return a.second > b.second;                             
    }
};

std::vector<int> AlgL(int,Graph&);
void randomGraph(int,int, Graph&);
void RMAT1(int,int, Graph&);

#endif  

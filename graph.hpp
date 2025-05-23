#ifndef GRAPH_HPP
#define GRAPH_HPP
#pragma once

#include <vector>
#include <utility>   
#include <string>

class Graph {
public:
    int n = 0;                                           
    int maxDist = 2e9;                                            
    std::vector<std::vector<std::pair<int,int>>> adj_lists;     

    void parse_graph(const std::string& filename);            
};

struct customCompare {
    bool operator()(const std::pair<int,int>& a,
                    const std::pair<int,int>& b) const noexcept {
        return a.second > b.second;                             
    }
};


#endif  
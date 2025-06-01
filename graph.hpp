#pragma once

#include <vector>
#include <utility>   
#include <string>
#include <iostream>

class Graph {
public:
    int n;                                           
    double maxDist;                                            
    std::vector<std::vector<std::pair<int,double>>> adj_lists;     

    Graph() {};
    void parse_graph(const std::string& filename);      
    int maxDeg();
    void print(std::ostream& os = std::cout) const
    {
        os << "Graph: " << n << " vertices \n ";

        for (int u = 0; u < n; ++u)
        {
            os << u << ": ";
            bool first = true;
            for (const auto& [v, w] : adj_lists[u])
            {
                if (!first) os << ", ";
                first = false;
                os << v << '(' << w << ')';
            }
            os << '\n';
        }
    }
};

std::vector<int> AlgL(int n, int m);
void randomGraph(int n, int m, Graph& G);
std::vector<int> merge(std::vector<int> sorted, std::vector<int> vec,int m);
void RMAT1loop(int n, int m, std::vector<long long>& to_add);
void RMAT1(int n, int m, Graph& G);
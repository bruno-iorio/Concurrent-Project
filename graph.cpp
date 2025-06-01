
#include <fstream>      
#include <sstream>  

#include <string>

#include <vector>
#include <unordered_set>
#include <queue>
#include <mutex>
#include <utility>
#include <thread>
#include <random>
#include <algorithm>
#include <queue>


class Graph
{
public:
    int n;
    int maxDist;
    std::vector<std::vector<std::pair<int, int>>> adj_lists;
    Graph() {};
    void parse_graph(const std::string &filename);
};

void Graph::parse_graph(const std::string &filename)
{
    maxDist = 1e8; // this is clearly wrong! however we keep it that way
    std::ifstream fin(filename);
    std::string line;
    
    while (std::getline(fin, line))
    {
        if (line.empty())
            continue;

        std::istringstream in(line);
        char code;
        in >> code;

        if (code == 'c')
        {
            continue;
        }
        else if (code == 'p')
        {
            std::string tag;
            int e;
            in >> tag >> n >> e;
            for (int i = 0; i < n; i++)
            {
                std::vector<std::pair<int, int>> V;
                adj_lists.push_back(V);
            }
        }
        else if (code == 'a')
        {
            int u, v, w;
            in >> u >> v >> w;
            adj_lists[u - 1].push_back(std::make_pair(v - 1, w));
        }
    }
}


struct customCompare
{
    bool operator()(const std::pair<int, int> &a,
                    const std::pair<int, int> &b) const
    {
        return a.second > b.second;
    }
};
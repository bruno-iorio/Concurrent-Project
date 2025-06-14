    #include "graph.hpp"
    #include <vector>
    #include <unordered_set>
    #include <queue>
    #include <mutex>
    #include <utility>
    #include <thread>
    #include <random>
    #include <algorithm>
    #include <queue>

    struct customCompare {
        bool operator()(const std::pair<int,double>& a,
                        const std::pair<int,double>& b) const noexcept {
            return a.second > b.second;                             
        }
    };

    class Dijkstra 
    {
    public:
        Graph graph;

        std::vector<double> tent;
        std::vector<bool> visited;

        
        Dijkstra(Graph graph)
        {
            this->graph = graph;
            this->tent.resize(graph.n);
            this->visited.resize(graph.n);
            std::fill(tent.begin(), tent.end(), 2e9);
            std::fill(visited.begin(), visited.end(), false);
        }

        void findShortest(int source)
        {
            std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, customCompare> pq;
            pq.push(std::make_pair(source, 0));

            while (!pq.empty())
            {
                std::pair<int, double> v = pq.top();
                pq.pop();

                // already optimised
                if (visited[v.first])
                    continue;

                // not yet
                tent[v.first] = std::min(tent[v.first], v.second);
                visited[v.first] = true;

                // relax
                for (auto u : graph.adj_lists[v.first])
                {
                    if (tent[u.first] > tent[v.first] + u.second)
                    {
                        pq.push(std::make_pair(u.first, tent[v.first] + u.second));
                    }
                }
            }
        }
    };

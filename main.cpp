#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <unordered_set>
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
            // initialise n vectors of outgoing edges
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

class DeltaStepping
{
public:
    int delta;
    Graph graph;
    std::vector<double> tent;
    std::vector<std::unordered_set<int>> buckets;

    DeltaStepping(Graph graph, double delta)
    {
        this->graph = graph;
        this->delta = delta;
        int n_buckets = graph.maxDist / delta + 1;
        this->buckets = std::vector<std::unordered_set<int>>(n_buckets);
        this->tent.resize(graph.n);
        std::fill(tent.begin(), tent.end(), 2e9);
    }

    void findShortest(int source)
    {
        relax(source, 0);
        int i;
        int maxiter = 10;
        while (!BucketsEmpty(i))
        {
            std::unordered_set<int> R;
            std::vector<std::pair<int, int>> Req;
            while (!buckets[i].empty())
            {
                Req = findRequests(buckets[i], true);
                for (auto vertex : buckets[i])
                    R.insert(vertex);
                buckets[i].clear();
                relaxRequests(Req);
            }
            Req = findRequests(R, false);
            relaxRequests(Req);
        }
    }

    // this works
    bool BucketsEmpty(int &non_empty_index)
    {
        for (int i = 0; i != buckets.size(); i++)
        {
            if (!buckets[i].empty())
            {
                non_empty_index = i;
                return false;
            }
        }
        return true;
    }

    // kind values:
    // 0 - heavy edges
    // 1 - light edges
    std::vector<std::pair<int, int>> findRequests(std::unordered_set<int> vertices, bool kind)
    {
        std::vector<std::pair<int, int>> reqs;

        for (auto vertex : vertices)
        {
            for (auto connected_vertex : graph.adj_lists[vertex])
            {
                if (kind && connected_vertex.second <= delta)
                    reqs.push_back(std::make_pair(connected_vertex.first, tent[vertex] + connected_vertex.second));
                else if (!kind && connected_vertex.second > delta)
                    reqs.push_back(std::make_pair(connected_vertex.first, tent[vertex] + connected_vertex.second));
            }
        }
        return reqs;
    }

    void relaxRequests(std::vector<std::pair<int, int>> reqs)
    {
        for (auto req : reqs)
        {
            relax(req.first, req.second);
        }
    }

    void relax(int vertex, int distance)
    {
        if (distance < tent[vertex])
        {
            // remove only if we are sure it is in a bucket
            if (tent[vertex] != 2e9)
                buckets[tent[vertex] / delta].erase(vertex);
            buckets[distance / delta].insert(vertex);
            tent[vertex] = distance;
        }
    }
};

struct customCompare
{
    bool operator()(const std::pair<int, int> &a,
                    const std::pair<int, int> &b) const
    {
        return a.second > b.second;
    }
};

class Dijkstra
{
public:
    Graph graph;
    std::vector<int> tent;
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
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, customCompare> pq;
        pq.push(std::make_pair(source, 0));

        while (!pq.empty())
        {
            std::pair<int, int> v = pq.top();
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

int main()
{
    Graph G;
    G.parse_graph("graphs/rome.gr");

    DeltaStepping alg1(G, 2);
    Dijkstra alg2(G);

    alg1.findShortest(0);
    alg2.findShortest(0);

    // do they find the same distance
    for (int i = 0; i < G.n; i++)
        std::cout << alg1.tent[i] - alg2.tent[i] << std::endl;
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <unordered_set>
#include <queue>
#include <random>
#include <algorithm>
#include <thread>
#include <chrono>


class Graph
{
public:
    int n;
    int maxDist;
    std::vector<std::vector<std::pair<int, double>>> adj_lists;
    Graph() {};
    void parse_graph(const std::string &filename);
    int maxDeg();
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
                std::vector<std::pair<int, double>> V;
                adj_lists.push_back(V);
            }
        }
        else if (code == 'a')
        {
            int u, v;
            double w;
            in >> u >> v >> w;
            adj_lists[u - 1].push_back(std::make_pair(v - 1, w));
        }
    }
}

int Graph::maxDeg(){
    int max_deg = 0;
    for(int i = 0; i < n; i++)
        max_deg = std::max(max_deg, int(adj_lists[i].size()));
    return max_deg;
}

//https://en.wikipedia.org/wiki/Reservoir_sampling
//algorithm L for sampling m values with no repetition from 1 to n.
// Sample m distinct integers from [1..n], returned in res (size m)
std::vector<int> AlgL(int n, int m) {
    std::random_device rd;
    std::mt19937_64 gen(rd());

    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::uniform_int_distribution<int> R(0, m-1);

    std::vector<int> res(m);

    double W = std::exp(std::log(U(gen)) / m);

    int i = m;
    while (true) {
        double u = U(gen);
        int s = int(std::floor(std::log(u) / std::log(1.0 - W))) + 1;
        i += s;
        if (i > n) break;

        int idx = R(gen);
        res[idx] = i;

        W *= std::exp(std::log(U(gen)) / m);
    }

    return res;
}

//https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
//this is how we get edge from a random index
//weights are U[0,1]
Graph randomGraph(int n, int m){
    Graph G;
    G.n = n;
    G.maxDist = 1e8;
    
    for (int i = 0; i < n; ++i) {
        std::vector<std::pair<int, double>> V;
        G.adj_lists.push_back(V);
    }
    std::vector<int> edges = AlgL(n * (n - 1) / 2 , m);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> U(0.0, 1.0);
    
    for(auto e : edges){
        int i = (1 + std::sqrt(1 + 8*e))/2;
        int j = e - i*(i-1)/2;
        std::pair<int,double> p1 = {j,U(gen)};
        std::pair<int,double> p2 = {i,U(gen)};
        std::cerr << e << '\n';
        std::cerr << i << " " <<  G.adj_lists.size() << '\n';
        std::cerr << j << " " << G.adj_lists.size() << '\n';
        G.adj_lists[i].push_back(p1);
        G.adj_lists[j].push_back(p2);
    }

    return G;
}

/*
Graph RMAT1{






};



Graph Grid{





};

*/






class DeltaSteppingSequential
{
public:
    double delta;
    Graph graph;
    std::vector<double> tent;
    std::vector<std::unordered_set<int>> buckets;

    DeltaSteppingSequential(Graph graph, double delta)
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
            std::vector<std::pair<int, double>> Req;
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
    std::vector<std::pair<int, double>> findRequests(std::unordered_set<int> vertices, bool kind)
    {
        std::vector<std::pair<int, double>> reqs;

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

    void relaxRequests(std::vector<std::pair<int, double>> reqs)
    {
        for (auto req : reqs)
        {
            relax(req.first, req.second);
        }
    }

    void relax(int vertex, double distance)
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
    bool operator()(const std::pair<int, double> &a,
                    const std::pair<int, double> &b) const
    {
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

class DeltaSteppingParallelStatic
{
public:
    double delta;
    int num_threads;
    Graph graph;
    std::vector<double> tent;
    std::vector<std::vector<std::unordered_set<int>>> buckets; //bucket - thread - set of elements
    std::vector<int> owner; 
    std::vector<std::vector<std::pair<int, double>>> NeighborsLight, NeighborsHeavy; //element - list of neighbors 
    std::vector<std::vector<std::vector<std::pair<int,double>>>> ReqLight, ReqHeavy; //thread - list of requests for this thread

    DeltaSteppingParallelStatic(const Graph& graph, double delta, int num_threads): 
        graph(graph),
        delta(delta),
        num_threads(num_threads),
        tent(graph.n, 2e9),
        buckets(
          graph.maxDist / delta + 1,
          std::vector<std::unordered_set<int>>(num_threads)
        ),
        owner(graph.n, 0),
        NeighborsLight(graph.n),
        NeighborsHeavy(graph.n),
        ReqLight(num_threads, std::vector<std::vector<std::pair<int, double>>>(num_threads)),
        ReqHeavy(num_threads, std::vector<std::vector<std::pair<int, double>>>(num_threads))
        {}

    //fills in the neighbours graph
    void fillNeighbors(){
        for(int i = 0; i < graph.n; i++){
            for(auto edge : graph.adj_lists[i]){
                if(edge.second <= delta)
                    NeighborsLight[i].push_back(edge);
                else
                    NeighborsHeavy[i].push_back(edge);
            }
        }
    }

    //assign owners to vertices
    void assignThreads(){
        int iter = 0;
        for (int t = 0; t < num_threads; t++) {
            int count = graph.n / num_threads + (t < (graph.n % num_threads) ? 1 : 0);
            for (int i = 0; i < count; ++i) {
                owner[iter++] = t;
            }
        }
    
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(owner.begin(), owner.end(), gen);
    }
    
    
    //thread, index of bucket
    void loop1(int t, int i){
        //generate Req_L(B_i_t) and Req_H(B_i_t)
        for(auto from : buckets[i][t]){
            //std::cout << "generating requests from " << from << std::endl;
            generateRequests(from, true);
            generateRequests(from, false);
        }
        buckets[i][t].clear();
    }

    void loop2(int t){
        Relax(true, t);
        for(int src = 0; src < num_threads; ++src) {
            ReqLight[src][t].clear();
        }
    }

    void loop3(int t){
        
        Relax(false, t);
        for(int src = 0; src < num_threads; ++src) {
            ReqHeavy[src][t].clear();
        }
    }

    //kind values:
    //0 - heavy requests
    //1 - light requests
    void generateRequests(int from, bool kind){
        if(kind){
            for(auto vertex : NeighborsLight[from]){
                ReqLight[owner[from]][owner[vertex.first]].push_back(std::make_pair(vertex.first, tent[from] + vertex.second));
            }
        }
        else{
            for(auto vertex : NeighborsHeavy[from]){
                ReqHeavy[owner[from]][owner[vertex.first]].push_back(std::make_pair(vertex.first, tent[from] + vertex.second));
            }

        }
    }

    //slightly different - relaxes the whole sublist for thread t instead of splitting it.
    void Relax(bool kind, int to_thread){
        if(kind){
            for(int from_thread = 0;  from_thread < num_threads; from_thread++){
                for(auto request : ReqLight[from_thread][to_thread]){
                    int vertex = request.first;
                    double propDist = request.second;
                    if(propDist < tent[vertex]){
                        tent[vertex] = propDist; 
                        int i = int(tent[vertex] / delta);
                        int j = int(propDist / delta);
                        if(tent[vertex] != 2e9)
                            buckets[i][to_thread].erase(vertex);
                        buckets[j][to_thread].insert(vertex);
                    }
                }
            }
        }
        else{
            for(int from_thread = 0;  from_thread < num_threads; from_thread++){
                for(auto request : ReqHeavy[from_thread][to_thread]){
                    int vertex = request.first;
                    double propDist = request.second;
                    if(propDist < tent[vertex]){
                        tent[vertex] = propDist; 
                        double i = tent[vertex] / delta;
                        double j = propDist / delta;
                        if(tent[vertex] != 2e9)
                            buckets[i][to_thread].erase(vertex);
                        buckets[j][to_thread].insert(vertex);
                    }
                }
            }
        }
    }

    bool BucketsEmpty(int &non_empty_index){
        for(int i = 0; i != buckets.size(); i++){
            for(int j = 0; j != buckets[i].size(); j++){
                if(!buckets[i][j].empty()){
                    non_empty_index = i;
                    return false;
                }
            }
        }
        return true;
    }

    bool BucketEmpty(int index){
        for(int j = 0; j != buckets[index].size(); j++){
            if(!buckets[index][j].empty()){
                return false;
            }
        }
        return true;
    }

    void findShortest(int source){
        //setup
        fillNeighbors();
        /*
        std::cout << "graph light edges " << std::endl;
        for(int i = 0; i < graph.n; i++){
            std::cout << i << ": ";
            for(auto j : NeighborsLight[i]){
                std::cout << j.first << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "graph heavy edges " << std::endl;
        for(int i = 0; i < graph.n; i++){
            std::cout << i << ": ";
            for(auto j : NeighborsHeavy[i]){
                std::cout << j.first << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        */
        assignThreads();

        /*
        std::cout << "owner threads" << std::endl;
        for(auto th : owner){
            std::cout << th << " ";
        }
        std::cout << std::endl;
        */

        //insert source to bucket 0
        int own0 = owner[source];
        buckets[0][own0].insert(source);

        //set tentative distance of source to 0
        tent[0] = 0;

        //start loops
        int index;
        std::vector<std::thread> workers(num_threads);
        while(!BucketsEmpty(index)){
            //std::cout << "entering nonempty bucket " << index << std::endl;
            //std::vector<std::thread> workers(num_threads);

            while(!BucketEmpty(index)){

                //std::cout << "entering loop 1" << std::endl;

                //1st loop
                for(int t = 0; t < num_threads; t++)
                    workers[t] = std::thread(&DeltaSteppingParallelStatic::loop1, this, t, index);
                for(int t = 0; t < num_threads; t++)
                    workers[t].join();

                //std::cout << "entering loop2" << std::endl;
                        
                //2nd loop
                for(int t = 0; t < num_threads; t++)
                    workers[t] = std::thread(&DeltaSteppingParallelStatic::loop2, this, t);
                for(int t = 0; t < num_threads; t++)
                    workers[t].join();
            }

            //std::cout << "entering loop 3" << std::endl;

            //3rd loop
            for(int t = 0; t < num_threads; t++)
                workers[t] = std::thread(&DeltaSteppingParallelStatic::loop3, this,  t);
            for(int t = 0; t < num_threads; t++)
                workers[t].join();
        }
    }
};


int main(int argc, char** argv) {
    /*
    std::vector<int> x = AlgL(20, 10);
    for(auto xs : x)
        std::cout << xs << " ";
    if(argc < 4) {
        std::cerr << "Usage: " << argv[0] << " threads delta option\n";
        return 1;
    }
    int threads = std::stoi(argv[1]);
    int delta   = std::stoi(argv[2]);
    int option  = std::stoi(argv[3]);

    Graph G;
    G.parse_graph("graphs/USA-road-d.BAY.gr");

    if(option == 0){
        Dijkstra alg(G);
        auto start = std::chrono::steady_clock::now();
        alg.findShortest(0);
        auto finish = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
        std::cout << elapsed << '\n';
    }
    else if(option == 1){
        DeltaSteppingSequential alg(G, delta);
        auto start = std::chrono::steady_clock::now();
        alg.findShortest(0);
        auto finish = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
        std::cout << elapsed << '\n';
    }
    else{
        DeltaSteppingParallelStatic alg(G, delta, threads);
        auto start = std::chrono::steady_clock::now();
        alg.findShortest(0);
        auto finish = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
        std::cout << elapsed << '\n';
        }

    
    */
    Graph G = randomGraph(10, 20);
    for(int i = 0; i < G.n; i++){
        std::cerr << i << ": ";
        for(auto v : G.adj_lists[i]){
            std::cerr << "(" << v.first << "," << v.second << ")" << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}

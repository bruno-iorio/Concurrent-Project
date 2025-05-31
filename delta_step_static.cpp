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

class DeltaSteppingSequential
{
public:
    double delta;
    Graph graph;
    std::vector<int> tent;
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
                buckets[(int) (tent[vertex] / delta)].erase(vertex);
            buckets[(int) (distance / delta)].insert(vertex);
            tent[vertex] = distance;
        }
    }
};



class DeltaSteppingParallelStatic
{
public:
    int delta;
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
          (int)(graph.maxDist / delta) + 1,
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
        while(!BucketsEmpty(index)){
            //std::cout << "entering nonempty bucket " << index << std::endl;
            std::vector<std::thread> workers(num_threads);

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


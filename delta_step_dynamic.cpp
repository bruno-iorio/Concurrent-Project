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


class DeltaSteppingSequentialDynamic {
public:

    DeltaSteppingSequentialDynamic(const Graph& g, int d)
        : delta(d), graph(g)
    {

        int nBuckets = graph.maxDist / delta + 1;
        buckets.resize(nBuckets);

        tent.assign(graph.n, 2e9);         
    }

    Graph graph;
    int delta;
    std::vector<int> tent;
    std::vector<std::unordered_set<int>> buckets; // here this grows on demand
    std::priority_queue<int,
                        std::vector<int>,
                        std::greater<int>> active; // min-heap of non-empty idx

    void ensureBucket(int idx) {
        if (idx >= static_cast<int>(buckets.size()))
            buckets.resize(idx + 1);
    }

    bool nextBucket(int& idx) {                   // O(log k)
        while (!active.empty()) {
            idx = active.top();
            if (!buckets[idx].empty()) 
                return false;
            active.pop();                         
        }
        return true;                             
    }


    void relax(int v, int d) {

        if (d >= tent[v]) 
            return;             

        if (tent[v] != 2e9) 
            buckets[tent[v] / delta].erase(v);

        tent[v] = d;
        int idx = tent[v] / delta;
        ensureBucket(idx);
        if (buckets[idx].empty()) 
            active.push(idx); 
        buckets[idx].insert(v);

    }

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

    void relaxRequests(const std::vector<std::pair<int,int>>& reqs) {
        for (auto [v,d] : reqs) 
            relax(v,d);
    }


    void findShortest(int s) {

        relax(s,0);
        int i;

        while (!nextBucket(i)) {

            std::unordered_set<int>  R;
            std::vector<std::pair<int,int>> Req;

            while (!buckets[i].empty()) {
                Req = findRequests(buckets[i], true);
                R.insert(buckets[i].begin(), buckets[i].end());
                buckets[i].clear();                
                relaxRequests(Req);
            }

            Req = findRequests(R, false);
            relaxRequests(Req);

        }
    }
};





class DeltaSteppingParallelDynamic
{
public:

    int delta;
    int num_threads;
    Graph graph;
    std::vector<int> tent;                                    
    std::vector<std::vector<std::unordered_set<int>>> buckets;  
    std::vector<int> owner;
    std::vector<std::vector<std::pair<int,double>>> neiLight , neiHeavy;
    std::vector<std::vector<std::vector<std::pair<int,double>>>> reqL, reqH;

    std::mutex resize_mtx;   
    std::mutex heap_mtx;     // protects globalActive
    std::priority_queue<int, std::vector<int>, std::greater<int>> globalActive;        

    DeltaSteppingParallelDynamic(const Graph& g, double d, int nThreads)
        : 
        graph(g),
        delta(d),
        num_threads(nThreads),
        tent(g.n, 2e9),
        buckets(1, std::vector<std::unordered_set<int>>(nThreads)),
        owner(g.n, 0),
        neiLight(g.n), neiHeavy(g.n),
        reqL(nThreads, std::vector<std::vector<std::pair<int, double>>>(nThreads)),
        reqH(nThreads, std::vector<std::vector<std::pair<int, double>>>(nThreads))
    { }

    void splitNeighbors()
    {
        for (int v = 0; v < graph.n; ++v)
            for (auto e : graph.adj_lists[v])
                if(e.second <= delta)
                    neiLight[v].push_back(e);
                else
                    neiHeavy[v].push_back(e);
    }

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

    void ensureBucket(int idx)
    {
        if (idx < (int)buckets.size()) return;
        std::lock_guard<std::mutex> g(resize_mtx);
        if (idx >= (int)buckets.size())   // check again after lock!
            buckets.resize(idx + 1, std::vector<std::unordered_set<int>>(num_threads));
    }

    inline void moveToBucket(int threadId, int v, int newD, int oldD)
    {
        if (oldD != 2e9) {
            int oldIdx = (int)(oldD / delta);
            buckets[oldIdx][threadId].erase(v);
        }

        int newIdx = (int)(newD / delta);
        ensureBucket(newIdx);

        if (buckets[newIdx][threadId].empty()) {      
            std::lock_guard<std::mutex> h(heap_mtx);
            globalActive.push(newIdx);
        }
        buckets[newIdx][threadId].insert(v);
    }

    void genReq(int from, bool kind)
    {
        if(kind){
            for(auto vertex : neiLight[from]){
                reqL[owner[from]][owner[vertex.first]].push_back(std::make_pair(vertex.first, tent[from] + vertex.second));
            }
        }
        else{
            for(auto vertex : neiHeavy[from]){
                reqH[owner[from]][owner[vertex.first]].push_back(std::make_pair(vertex.first, tent[from] + vertex.second));
            }

        }
    }

    void loop1(int th, int idx)  
    {
        for (int v : buckets[idx][th]) 
        {
            genReq(v,true); 
            genReq(v,false); 
        }

        buckets[idx][th].clear();
    }

    void relaxReqBuf(bool light, int dstT)
    {
        auto& buf = light ? reqL : reqH;
        for (int srcT = 0; srcT < num_threads; ++srcT)
            for (auto [v,d] : buf[srcT][dstT]) {
                if (d < tent[v]) {
                    int old = tent[v];
                    tent[v]  = d;
                    moveToBucket(dstT, v, d, old);
                }
            }
        for (int srcT = 0; srcT < num_threads; ++srcT)
            buf[srcT][dstT].clear();
    }

    void loop2(int t) 
        {relaxReqBuf(true , t);}

    void loop3(int t) 
        {relaxReqBuf(false, t);}

    bool nextBucket(int &idx)
    {
        std::lock_guard<std::mutex> h(heap_mtx);
        while (!globalActive.empty()) {
            idx = globalActive.top();
            /* is the whole bucket[idx] empty across all threads? */
            bool empty = true;
            for (int t = 0; t < num_threads && empty; ++t)
                empty = empty && buckets[idx][t].empty();
            if (!empty) 
                return true;              
            globalActive.pop();                    // stale index, discard
        }
        return false;                              
    }

    void findShortest(int s)
    {
        splitNeighbors();
        assignThreads();

        tent[s] = 0;
        int t0  = owner[s];
        buckets[0][t0].insert(s);
        globalActive.push(0);

        std::vector<std::thread> pool(num_threads);
        int idx;

        while (nextBucket(idx)) {

            while (true) {                 // LIGHT phase loops 1+2
                std::vector<std::thread> pool(num_threads);

                // generate requests from current bucket
                for (int t = 0; t < num_threads; ++t)
                    pool[t] = std::thread(&DeltaSteppingParallelDynamic::loop1,
                                          this, t, idx);
                for (auto& th : pool) th.join();

                // relax light requests
                for (int t = 0; t < num_threads; ++t)
                    pool[t] = std::thread(&DeltaSteppingParallelDynamic::loop2,
                                          this, t);
                for (auto& th : pool) th.join();

                bool empty = true;
                for (int t = 0; t < num_threads && empty; ++t)
                    empty = empty && buckets[idx][t].empty();
                if (empty)
                    break;
            }

            // heavy requests
            for (int t = 0; t < num_threads; ++t)
                pool[t] = std::thread(&DeltaSteppingParallelDynamic::loop3,
                                      this, t);
            for (auto& th : pool) 
                th.join();
        }
    }
};
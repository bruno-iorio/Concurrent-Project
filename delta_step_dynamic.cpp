#include "graph.hpp"
#include <vector>
#include <queue>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <random>
#include <algorithm>
#include <atomic>
#include <barrier>
#include <iostream>


class DeltaSteppingSequentialDynamic {

    const int INF = 2e8;
    const Graph &G;
    int delta;
    const int delta_max;
    std::vector<int> tent;
    std::vector<std::unordered_set<int>> buckets;
    
    std::priority_queue<int,std::vector<int>,std::greater<int>> active;

    void ensureBucket(int idx) {
        if (idx >= buckets.size()) 
            buckets.resize(idx+1); 
    }

    bool nextBucket(int &idx) {
        while (!active.empty()) 
        {
            idx = active.top();
            if (!buckets[idx].empty()) 
                return true;
            active.pop();
        }
        return false;
    }

    void rebuild_buckets() {

        std::vector<int> live;
        for (std::size_t idx = 0; idx < buckets.size(); ++idx)
            for (int v : buckets[idx]) 
                live.push_back(v);

        buckets.assign(G.maxDist / delta + 1, {});

        while (!active.empty()) 
            active.pop();

        for (int v : live) {
            int idx = tent[v] / delta;
            if (buckets[idx].empty()) 
                active.push(idx);
            buckets[idx].insert(v);
        }
    }


    void relax(int v,int d) {
        if (d >= tent[v]) 
            return;
        if (tent[v] != INF) 
            buckets[tent[v]/delta].erase(v);
        tent[v] = d;
        int idx = d / delta;  

        ensureBucket(idx);

        if (buckets[idx].empty()) 
            active.push(idx);

        buckets[idx].insert(v);
    }

    std::vector<std::pair<int,int>> genReq(const std::unordered_set<int>& S,bool light) {
        std::vector<std::pair<int,int>> R;  
        R.reserve(S.size()*4);
        for (int u : S)
            for (auto [v,w] : G.adj_lists[u])
                if ((light && w<=delta) || (!light && w>delta))
                    R.emplace_back(v, tent[u]+w);
        return R;
    }

    void relaxReq(const std::vector<std::pair<int,int>>& R) {
        for (auto [v,d] : R) 
            relax(v,d);
    }

public:
    DeltaSteppingSequentialDynamic(const Graph& g,int d)
        : G(g), delta(d), delta_max(std::max(1, g.maxDist/8)), tent(g.n,INF), buckets(g.maxDist/delta + 1) {}

    const std::vector<int>& distances() {
        return tent; 
    }

    void findShortest(int s) {
        relax(s,0); 
        int idx;
        while (nextBucket(idx)) {
            std::unordered_set<int> R;
            int lightRounds = 0;
            std::size_t cntLight = 0, cntHeavy = 0;

            // light edges
            while (!buckets[idx].empty()) {
                ++lightRounds;
                auto reqL = genReq(buckets[idx], true);
                cntLight += reqL.size();
                R.insert(buckets[idx].begin(), buckets[idx].end());
                buckets[idx].clear();
                relaxReq(reqL);
            }

            // heavy edges
            auto reqH = genReq(R, false);
            cntHeavy = reqH.size();
            relaxReq(reqH);

            if (lightRounds > 20 && delta < delta_max) {
            delta *= 2;
            rebuild_buckets();
        }

        }
    }
};


class DeltaSteppingParallelDynamic {
    const int INF = 2e9;
    const int MAX_REBUILDS;
    const int delta_update;


    enum Phase {IDLE, GEN_REQ_LIGHT, RELAX_LIGHT, RELAX_HEAVY, EXIT};

    const Graph &G;
    int delta;
    const int T;             
    int rebuild_cnt = 0;         
    const int delta_max;          

    std::vector<int> tent;
    std::vector<char> inBucket;  // vertex live flag
    std::vector<int>  owner;  // vertex to thread

    std::vector<std::vector<std::vector<int>>> buckets;
    std::priority_queue<int,std::vector<int>,std::greater<int>> activeHeap;
    std::mutex heap_mtx, resize_mtx;

    // pre-split neighbor lists to skip weight check 
    std::vector<std::vector<std::pair<int,int>>> neiLight, neiHeavy;
    std::vector<std::vector<std::vector<std::pair<int,int>>>> reqL, reqH;

    std::vector<std::thread> pool;
    std::barrier<> phase_barrier;
    std::atomic<Phase> phase{IDLE};
    std::atomic<int> curBucket{0};

    void splitNeighbors() {
        neiLight.resize(G.n);  
        neiHeavy.resize(G.n);
        for (int u = 0; u < G.n; ++u)
            for (auto [v,w] : G.adj_lists[u])
        if (w <= delta) 
            neiLight[u].push_back({v, w});
        else 
            neiHeavy[u].push_back({v, w});
    }

    void assignThreads() {
        owner.resize(G.n);
        int base = G.n / T, extra = G.n % T, idx = 0;
        for (int t = 0; t < T; ++t) {
            int cnt = base + (t<extra);
            for (int i = 0; i < cnt; ++i) 
                owner[idx++] = t;
        }
        // shuffle threads to avoid imbalanced loads
        std::shuffle(owner.begin(), owner.end(), std::mt19937{std::random_device{}()});
    }
    void startWorkers() {
        pool.resize(T);
        for (int t = 0; t < T; ++t)
            pool[t] = std::thread(&DeltaSteppingParallelDynamic::worker, this, t);
    }


    void rebuild_all_for_new_delta(int newDelta)
{
    delta = newDelta;                                  

    for (int u = 0; u < G.n; ++u) {
        neiLight[u].clear(); 
        neiHeavy[u].clear();
        for (auto [v,w] : G.adj_lists[u])
            if (w <= delta) 
                neiLight[u].push_back({v, w});
            else 
                neiHeavy[u].push_back({v, w});
    }

    std::vector<int> live;
    for (auto &perIdx : buckets)
        for (auto &vecT : perIdx)
            live.insert(live.end(), vecT.begin(), vecT.end());

    buckets.assign(G.maxDist/delta + 1, std::vector<std::vector<int>>(T));
    while (!activeHeap.empty()) 
        activeHeap.pop();

    for (int v : live) {
        int tid = owner[v];
        int idx = tent[v] / delta;
        auto &vec = buckets[idx][tid];
        if (vec.empty()) { 
            std::lock_guard l(heap_mtx); 
            activeHeap.push(idx); 
        }
        vec.push_back(v);
    }
}

    void ensureBucket(int idx) {
        if (idx < (int)buckets.size()) 
            return;
        std::lock_guard lg(resize_mtx);
        if (idx >= (int)buckets.size())
            buckets.resize(idx+1,std::vector<std::vector<int>>(T));
    }

    void insertVertex(int tid, int v, int newD, int oldD) {
        if (oldD != INF) {
            int oldIdx = oldD / delta;
            auto &vec = buckets[oldIdx][tid];
            vec.erase(std::remove(vec.begin(),vec.end(),v), vec.end());
        }
        int idx = newD / delta;
        ensureBucket(idx);
        auto &vec = buckets[idx][tid];
        if (vec.empty()) { 
            std::lock_guard l(heap_mtx); 
            activeHeap.push(idx); 
        }
        vec.push_back(v); 
        inBucket[v] = 1;
    }
    bool bucketEmpty(int idx) const {
        for (int t = 0; t < T; ++t)
            if (!buckets[idx][t].empty()) 
                return false;
        return true;
    }

    bool nextBucket(int &idx) {
        std::lock_guard l(heap_mtx);
        while (!activeHeap.empty()) {
            idx = activeHeap.top();
            if (!bucketEmpty(idx)) 
                return true;     
            activeHeap.pop();                      
        }
        return false;
    }

    // generate requests from bucket
    void loop1(int tid, int idx) {       
        for (int u : buckets[idx][tid]) {
            for (auto [v,w] : neiLight[u])
                reqL[tid][owner[v]].push_back({v, tent[u]+w});
            for (auto [v,w] : neiHeavy[u])
                reqH[tid][owner[v]].push_back({v, tent[u]+w});
        }
        buckets[idx][tid].clear();
    }
    
    void relaxBuf(bool light, int dstT)
    {
        auto &buffer = (light ? reqL : reqH);

        for (int srcT = 0; srcT < T; ++srcT)
        {
            auto &reqList = buffer[srcT][dstT];

            for (const auto &req : reqList)
            {
                int v = req.first;   
                int d = req.second;  

                if (d < tent[v])                
                {
                    int oldDist = tent[v];
                    tent[v] = d;                  
                    insertVertex(dstT, v, d, oldDist);   
                }
            }

            reqList.clear();
        }
    }

    void loop2(int t) {
        relaxBuf(true, t);
    }

    void loop3(int t) {
        relaxBuf(false,t);
    }


    void startPhase(Phase p) {
        phase.store(p, std::memory_order_relaxed);  
        phase_barrier.arrive_and_wait();         
    }

    void worker(int tid) {
        bool run = true;
        while (run) {
            phase_barrier.arrive_and_wait();         
            switch (phase.load(std::memory_order_relaxed)) {
            case GEN_REQ_LIGHT: 
                loop1(tid, curBucket.load());          
                break;
            case RELAX_LIGHT: 
                loop2(tid);                            
                break;
            case RELAX_HEAVY:    
                loop3(tid);          
                break;
            case EXIT:           
                run = false;         
                break;
            }
            phase_barrier.arrive_and_wait();       
        }
    }

public:
    DeltaSteppingParallelDynamic(const Graph& g,int d,int threads, int light_threshold = 60, int max_rebuilds = 5)
        : G(g), delta(d), T(threads), delta_max(std::max(1, g.maxDist/8)),
          delta_update(light_threshold), MAX_REBUILDS(max_rebuilds),
          tent(g.n,INF), inBucket(g.n,0),
          buckets(1,std::vector<std::vector<int>>(threads)),
          reqL(threads,std::vector<std::vector<std::pair<int,int>>>(threads)),
          reqH(threads,std::vector<std::vector<std::pair<int,int>>>(threads)),
          phase_barrier(threads+1)                // +1 main
    {
        splitNeighbors();
        assignThreads();
        startWorkers();
    }

    ~DeltaSteppingParallelDynamic() {  
        for (auto &t : pool) 
            t.join();
    }

    const std::vector<int>& distances() {
        return tent;
    }

    void findShortest(int s) {
        tent[s] = 0;
        insertVertex(owner[s], s, 0, INF);

        int idx;
        while (nextBucket(idx)) {
            curBucket.store(idx,std::memory_order_relaxed);

            int lightRounds = 0;

            while (true) {
                ++lightRounds;                 
                startPhase(GEN_REQ_LIGHT);   
                phase_barrier.arrive_and_wait();
                startPhase(RELAX_LIGHT);   
                phase_barrier.arrive_and_wait();
                if (bucketEmpty(idx))       
                    break;                
            }

            startPhase(RELAX_HEAVY);
            phase_barrier.arrive_and_wait();

            // number of lightRounds we want before delta changes is dependant on graph size
            if (lightRounds > delta_update && delta < delta_max && rebuild_cnt < MAX_REBUILDS) {
                rebuild_all_for_new_delta(delta * 2);
                ++rebuild_cnt;
            }
        }

        startPhase(EXIT);    

    }
};

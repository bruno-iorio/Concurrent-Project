//we used chatgpt when building 
/* 
g++ -std=c++20 -O3 -pthread graph.cpp dijkstra.cpp delta_step_static.cpp delta_step_dynamic.cpp alg_compare.cpp -o alg_compare
*/
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>

#include "dijkstra.cpp"
#include "delta_step_static.cpp"   
#include "delta_step_dynamic.cpp"  

/* ------------------------------------------------------------------ */
/* helpers                                                            */
/* ------------------------------------------------------------------ */
using Clock = std::chrono::steady_clock;

template <typename Algo>
double timeAlgo(const char* name, Algo&& algo, int source)
{
    auto t0 = Clock::now();
    algo.findShortest(source);
    auto t1 = Clock::now();
    double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << std::left << std::setw(32) << name
              << " : " << std::fixed << std::setprecision(2)
              << ms << " ms\n";
    return ms;
}

template <typename DistA, typename DistB>
void checkEqual(const std::vector<DistA>& a,
                const std::vector<DistB>& b,
                const char* what)
{
    if (a.size() != b.size()) {
        std::cerr << "Size mismatch in " << what << "\n";
        std::abort();
    }
    for (size_t i = 0; i < a.size(); ++i)
        if (static_cast<long long>(a[i]) !=
            static_cast<long long>(b[i])) {
            std::cerr << "Mismatch at vertex " << i
                      << " (" << what << ")\n";
            std::abort();
        }
}

/* ------------------------------------------------------------------ */
/* main                                                               */
/* ------------------------------------------------------------------ */
int main(int argc, char* argv[])
{
    const std::string file =
        (argc > 1) ? argv[1] : "graphs/rome.gr";
    const int  SRC        = 0;
    const double  DELTA   = 0.5;                 // same δ everywhere
    const int  THREADS    = 10;

    /* ---------- load graph ------- --------------------------------- */
    Graph G;

    wideWeightRandomGraph(1<<17, (1<<17) * 100, 0, 1, G);

   /* try {
        G.parse_graph(file);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }  */

    std::cout << "\nBenchmark on " << file
              << "  (|V| = " << G.n
              << ")\n\n";

    /* ---------- construct algorithms ------------------------------ */
    Dijkstra                      a_dij   (G);
    //DeltaSteppingSequential       a_seqS  (G, DELTA);
    //DeltaSteppingSequentialDynamic a_seqD (G, DELTA);
    DeltaSteppingParallelStatic   a_parS  (G, DELTA, THREADS);
    DeltaSteppingParallelDynamic  a_parD  (G, DELTA, THREADS);

    /* ---------- run + time ---------------------------------------- */
    timeAlgo("Binary-heap Dijkstra",           a_dij,  SRC);
    //timeAlgo("Δ-Stepping sequential (static)", a_seqS, SRC);
    //timeAlgo("Δ-Stepping sequential (dynamic)",a_seqD, SRC);
    timeAlgo("Δ-Stepping parallel  (static)",  a_parS, SRC);
    timeAlgo("Δ-Stepping parallel  (dynamic)", a_parD, SRC);

    /* ---------- correctness gate ---------------------------------- */
    //checkEqual(a_dij.tent, a_seqS.tent, "Dijkstra vs seq-static");
 //   checkEqual(a_ dij.tent, a_seqD.tent, "Dijkstra vs seq-dynamic");
    checkEqual(a_dij.tent, a_parS.tent, "Dijkstra vs par-static");
    checkEqual(a_dij.tent, a_parD.distances(),  "Dijkstra vs par-dynamic");

    std::cout << "\nAll algorithms returned identical distances ✅\n";
    return 0;
}


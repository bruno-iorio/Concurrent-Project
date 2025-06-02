// alg_compare.cpp
// -------------------------------------------------------------
// Build example:
//
// g++ -std=c++17 -O3 -pthread \
//     graph.cpp dijkstra.cpp delta_stepping_static.cpp \
//     delta_stepping_dynamic.cpp alg_compare.cpp -o alg_compare
// -------------------------------------------------------------
//CHATPGT WAS USED WHEN BUILDING THIS FILE

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
#include "graph.hpp"

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
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
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


void benchmark_roadmap(const std::string filename, int option, int delta, int threads, int light = 60 , int rebuilds = 5){
	Graph G;
	G.parse_graph(filename);
	double time;
	if(option == 0){
    		Dijkstra alg(G);
    		time = timeAlgo("Dijkstra", alg, 0);
	}
	else if(option == 1){
    		DeltaSteppingSequential alg(G,delta);
    		time = timeAlgo("dss", alg, 0);
	}
	else if(option == 2){
    		DeltaSteppingParallelStatic alg(G,delta,threads);
    		time = timeAlgo("dsps", alg, 0);
	}
	else if(option == 3){
	DeltaSteppingParallelDynamic alg(G,delta,threads,light,rebuilds);
    		time = timeAlgo("dspd", alg, 0);
	}
	std::cout << time << "\n";
}

void benchmark_random_graph(int c, int option, double delta, int threads, int avg, int light=60, int rebuilds=5){
	Graph G;
	int n = 1 << 17;
	int m = n * c; // n times density
	wideWeightRandomGraph(n,m,0,10,G);
	double time = 0;
	for(int i = 0 ; i != avg; i++){
		if(option == 0){
			Dijkstra alg(G);
			time += timeAlgo("Dijkstra", alg, 0);
		}
		else if(option == 1){
			DeltaSteppingSequential alg(G,delta);
			time += timeAlgo("dsss", alg, 0);
		}
		else if(option == 2){
			DeltaSteppingParallelStatic alg(G,delta,threads);
			time += timeAlgo("dsps", alg, 0);
		}
		else if(option == 3){
			DeltaSteppingParallelDynamic alg(G,delta,threads,light,rebuilds);
			time += timeAlgo("dspd", alg, 0);
		}
	}
	std::cout << time / (double) avg << '\n';
}



int main(int argc, char* argv[])
{
    if(argc != 9) {
        std::cout << "Usage: " << argv[0] << " threads delta option average density light rebuilds road_map\n";
        return 1;
    }
    int          threads = std::stoi(argv[1]);
    double         delta = std::stod(argv[2]);
    int           option = std::stoi(argv[3]);
    int          average = std::stoi(argv[4]);
    int                c = std::stoi(argv[5]);
    int            light = std::stoi(argv[6]);
    int         rebuilds = std::stoi(argv[7]);
    std::string filename = argv[8];
    
    if(!filename.ends_with(".gr")) benchmark_random_graph(c,option,delta,threads,average,light,rebuilds);
    else benchmark_roadmap(filename,option,delta,threads,light,rebuilds);
    return 0;
}

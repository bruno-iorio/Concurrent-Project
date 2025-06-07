
//CHATPGT WAS USED WHEN BUILDING THIS FILE

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>
#include <sys/resource.h>


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

double get_cpu_time() {// this implementation might not have the same precision as steady_clock -> good to keep in mind
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        perror("getrusage failed");
        return 0.0;
    }
    double user = usage.ru_utime.tv_sec * 1000.0 + usage.ru_utime.tv_usec / 1000.0;
    double sys  = usage.ru_stime.tv_sec * 1000.0 + usage.ru_stime.tv_usec / 1000.0;
    return user + sys;
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


void benchmark_roadmap(const std::string filename, int option, double delta, int threads, int light = 60 , int rebuilds = 5){
	Graph G;
	G.parse_graph(filename);
	double time;
	double usage;
	
	double cpu_begin = get_cpu_time();
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
   	double cpu_end = get_cpu_time();
	
	usage = (cpu_end - cpu_begin) / time;  
	

	std::cout << time << " " << usage << "\n";
}

void benchmark_random_graph(int c, int option, double delta, int threads,  int light=60, int rebuilds=5){
	Graph G;
	int n = 1 << 17;
	int m = n * c; // n times density
	
	wideWeightRandomGraph(n,m,0,1,G);
	double time = 0;
	double usage;
	double cpu_begin = get_cpu_time();
	
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

   	double cpu_end = get_cpu_time();
	usage = (cpu_end - cpu_begin) / time;  
	std::cout << time  << " " << usage << '\n';
}

int main(int argc, char* argv[])
{
    if(argc != 8) {
        std::cout << "Usage: " << argv[0] << " threads delta option density light rebuilds road_map\n";
        return 1;
    }
    int          threads = std::stoi(argv[1]);
    double         delta = std::stod(argv[2]);
    int           option = std::stoi(argv[3]);
    int                c = std::stoi(argv[4]);
    int            light = std::stoi(argv[5]);
    int         rebuilds = std::stoi(argv[6]);
    std::string filename = argv[7];
    
    if(!filename.ends_with(".gr")) benchmark_random_graph(c,option,delta,threads,light,rebuilds);
    else benchmark_roadmap(filename,option,delta,threads,light,rebuilds);
    return 0;
}

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
#include <cassert>
#include "graph.hpp"
#include "delta_step_dynamic.cpp"
#include "delta_step_static.cpp"
#include "dijkstra.cpp"

int main(int argc, char** argv) {
    if(argc != 9) {
        std::cout << "Usage: " << argv[0] << " threads delta option option_graph average density light rebuilds\n";
        return 1;
    }
    int      threads = std::stoi(argv[1]);
    double   delta   = std::stod(argv[2]);
    int      option  = std::stoi(argv[3]);
    int option_graph = std::stoi(argv[4]);
    int      average = std::stoi(argv[5]);
    int            c = std::stoi(argv[6]);
    int        light = std::stoi(argv[6]);
    int     rebuilds = std::stoi(argv[6]);
    
    int n = 1 << 3;
    int m = c * n;

    std::vector<Graph> G(average);
    std::vector<std::thread> workers(average);
    
    for (int i = 0; i !=average; i++){
    	randomGraph(n,m,G[i]);
    }

    if (option_graph == 0){
      for (int i = 0; i != average; i++){
        workers[i] = std::thread(randomGraph,n,m,std::ref(G[i]));
      }
    }

    if (option_graph == 1){
      for (int i = 0; i != average; i++){
        workers[i] = std::thread(RMAT1,n,m,std::ref(G[i]));
      }
    }
    for(int i = 0; i != average; i++){
      workers[i].join();
    }
    int64_t total_time = 0;
    for(int i=0; i!= average; i++){
      if(option == 0){
          Dijkstra alg(G[i]);
          auto start = std::chrono::steady_clock::now();
          alg.findShortest(0);
          auto finish = std::chrono::steady_clock::now();
          auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
          total_time += elapsed;

      }

      else if(option == 1){
          DeltaSteppingSequentialDynamic alg(G[i], delta);
          auto start = std::chrono::steady_clock::now();
          alg.findShortest(0);
          auto finish = std::chrono::steady_clock::now();
          auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
          total_time += elapsed;
      }
      
      else if (option == 2){
          DeltaSteppingSequential alg(G[i], delta);
          auto start = std::chrono::steady_clock::now();
          alg.findShortest(0);
          auto finish = std::chrono::steady_clock::now();
          auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
          total_time += elapsed;
      }
      
      else if(option == 3){
          DeltaSteppingParallelDynamic alg(G[i], (int) delta, threads);
          auto start = std::chrono::steady_clock::now();
          alg.findShortest(0);
          auto finish = std::chrono::steady_clock::now();
          auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
          total_time += elapsed;
      }
      
      else if (option == 4){
          DeltaSteppingParallelStatic alg(G[i], delta, threads);
          auto start = std::chrono::steady_clock::now();
          alg.findShortest(0);
          auto finish = std::chrono::steady_clock::now();
          auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
          total_time += elapsed;
      }
    }
    std::cout <<total_time<< '\n';
    return 0;
}

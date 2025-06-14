
#pragma once

#include <fstream>      
#include <sstream>  
#include <string>
#include <vector>
#include <unordered_set>
#include <queue>
#include <mutex>
#include <utility>
#include <thread>
#include <random>
#include <algorithm>
#include <queue>
#include <set>

#include "graph.hpp"


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
            for (int i = 0; i < n; i++)
            {
                std::vector<std::pair<int, double>> V;
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

int Graph::maxDeg()
{
    int max_deg = 0;
    for (int i = 0; i < n; i++)
        max_deg = std::max(max_deg, int(adj_lists[i].size()));
    return max_deg;
}

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
        int s = (std::floor(std::log(u) / std::log(1.0 - W))) + 1;
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

void randomGraph(int n, int m, Graph& G){
    G.n = n;
    for (int i = 0; i < n; ++i) {
        std::vector<std::pair<int, double>> V;
        G.adj_lists.push_back(V);
    }
    std::vector<int> edges = AlgL(n * (n - 1) / 2 , m);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for(auto e : edges){
        int i = std::floor((1 + std::sqrt(1 + 8*e)) / 2) ;
        int j = e - i*(i-1)/2;
        i--;
        std::pair<int,double> p1 = std::make_pair(j,U(gen));
        std::pair<int,double> p2 = std::make_pair(i,U(gen));
        G.adj_lists[i].push_back(p1);
        G.adj_lists[j].push_back(p2);
        G.maxDist += p1.second + p2.second;
    }
}

std::vector<int> merge(std::vector<int> sorted, std::vector<int> vec,int m){
  std::vector<int> res(sorted.size() + vec.size());
  std::sort(vec.begin(),vec.end());
  int i1= 0, i2= 0;
  while(i1 != sorted.size() && i2 != vec.size() && res.size() < m){
  	if(res.size() == 0){
	if (sorted[i1] < vec[i2]){
		res.push_back(sorted[i1]);
		i1++;
	}
	else if(sorted[i1] > vec[i2]){
		res.push_back(vec[i2]);
		i2++;
	}
	else{
		res.push_back(vec[i2]);
		i2++;
		i1++;
	}
     }
     else{
	if (sorted[i1] != *(res.end() - 1) && sorted[i1] < vec[i2]){
		res.push_back(sorted[i1]);
		i1++;
	}
	else if( vec[i2] != *(res.end() - 1) && sorted[i1] > vec[i2]){
		res.push_back(vec[i2]);
		i2++;
	}
	else if(sorted[i1] == *(res.end() -1)){
		i1++;
	}
	else if(vec[i1] == *(res.end()-1)){
		i2++;
	}
	else if(sorted[i1] == vec[i2]){
		res.push_back(vec[i2]);
		i2++;
		i1++;
	}
     }
  } 
  return res;
}

void RMAT1loop(int n, int m, std::vector<long long>& to_add){
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> U(0.0, 1.0);
  double pa = 0.45;
  double pb = 0.22;
  double pc = 0.22;
  int k = 0;
  while(k < m){
    int col = 0, row = 0;
    int size = n / 2;
    while (size > 0){
      double p = U(gen);
      if (p < pa){ //upper-left quadrant
      	  size = size / 2;
          continue;
      }
      else if (p < pa + pb){
        col += size;
      }
      else if (p < pa + pb + pc){
        row += size;
      }
      else{
        col += size;
        row += size;
      }
      size = size / 2;
    }
    long long i = std::max(col, row);
    long long j = std::min(col, row);
    long long e = i * n + j;
    if (i!= j){
    	to_add.push_back(e);
    	k++;
    }
    }
}

void RMAT1(int n, int m, Graph& G){
  int n_threads = 10;
  G.n = n;
  G.maxDist = 1e8;
  int elements_per_thread = m / n_threads;
  std::vector<std::thread> workers(n_threads);
  std::vector<std::vector<long long>> to_add(n_threads);
  std::unordered_set<long long> edges;
  
  for (int i = 0; i < n; ++i) {
    std::vector<std::pair<int, double>> V;
    G.adj_lists.push_back(V);
  }
  
  while(edges.size() < m){
  for (int i = 0; i != n_threads; i++){
  	workers[i] = std::thread(RMAT1loop, n, elements_per_thread,std::ref(to_add[i]));
  }
  for (int i = 0; i != n_threads; i++){
  	workers[i].join();
  }
  for(int i = 0 ; i!= n_threads ; i++){
  	edges.insert(to_add[i].begin(),to_add[i].end());
	to_add[i].clear();
	if(edges.size() >= m) break; 
  }
}

std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_real_distribution<double> U(0.0, 1.0);
int count = 0;
for(auto it : edges){
	if(count >= m) break;
	count++;
        int j = it % n ;
        int i = it / n ;
	G.adj_lists[i].push_back(std::make_pair(j, U(gen)));
	G.adj_lists[j].push_back(std::make_pair(i, U(gen)));
}
}

// credit: CHATGPT
void wideWeightRandomGraph(int n, int m, double min_w, double max_w, Graph& G) {
    G.n = n;
    G.adj_lists.assign(n, {});
    G.maxDist = 0;
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<int> node_dist(0, n-1);
    std::uniform_real_distribution<double> weight_dist(min_w, max_w);
    std::set<std::pair<int, double>> edges;
    while (edges.size() < (size_t)m) {
        int u = node_dist(gen), v = node_dist(gen);
        if (u == v) continue;
        int a = std::min(u, v), b = std::max(u, v);
        if (!edges.insert({a, b}).second) continue;
        double w = weight_dist(gen);
        G.adj_lists[a].emplace_back(b, w);
        G.adj_lists[b].emplace_back(a, w);
        G.maxDist += w;
    }
}
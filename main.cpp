#include <iostream>
#include <vector>
#include <limits>

void initializeGraph(std::vector<std::vector<int>>& Edges,int& n, int& maxEdge){
  for(int i = 0; i < n; i++){
    for(int j = 0; j != i; j++){
      std::cout << "Input the weight of Edge (-1 means unreachable)" << i << " and " << j " ";
      std::cin >> Edges[i][j];
      Edges[j][i] = Edges[i][j];
      if (Edges[i][j] > maxEdge){
        maxEdge = Edges[i][j]; 
      }
    }
  }
}


class DeltaStepping{ // sequential version of the algorithm
  public:
    DeltaStepping(std::vector<std::vector<int>> graph, int maxEdge, double delta){
      this->n = graph.size();
      this->maxEdge = maxEdge;
      this->delta = delta;

      int n_buckets = maxEdge / delta  + 1;
      B = std::vector<std::vector<int>>(n_buckets);

      tent = std::vector<int>(n,-1);

    }
    void findShortest(int source, int destination){
      relax(source, 0);
      std::vector<int> light, heavy;
      int i;
      while (!BucketEmpty(i)){
        std::vector<int> R;
        std::vector<std::pair<int,int>> Req;
        while (!Buckets[i].empty()){
          Req  = findRequests(B[i],light);
          for(int e : B[i]){
            if (!std::find(B[i].begin(),B[i].end(),e)) R.push_back(e); 
          }
            B[i].clear();
            relaxRequest(Req);
          }
        Req = findRequests(R,heavy);
        relaxRequest(Req);
        }
      }
    }

    std::vector<std::pair<int,int>> findRequests(std::vector<int> V, std::vector<int> kind){
      res = std::vector<std::pair<int,int>>;
      for(int v : V){
        for (int w = 0; w != n; w++){
          if(graph[w][v] >= 0){
            res.push_back({w, tent[v] + graph[w][v]});
          }
        }
      }
      return res;
    }
    
    void relaxRequest(std::vector<std::pair<int,int>> Req){
      for(auto request : Req){
        relax(request.first,request.second);
      }
    }

    bool BucketEmpty(int& non_empty_index){
      for (int i = 0; i != B.size(); i++){
        if (B[i].empty()) {
          non_empty_index = i;
          return false;
        }
      }
      return true;
    }

    void relax(int w, int x){
      if (x < tent[w]){
        B[tent(w)/delta].erase(B[tent(w) / delta].begin(), B[tent(w)/delta].end(), w);
        B[x / delta].push_back(w);
        tent[w] = x;
      }
    }
    int delta; 
    std::vector<double> tent;
    std::vector<std::vector<int>> B;
}

void Dijkstra(std::vector<std::vector<int>> Graph, int source, int destination); // to be implemented

int main(){
  int n; // number of nodes
  std::cout << " Input the number of Edges: ";
  std::cin >> n;
  std::vector<std::vector<int>> Graph(n, std::vector<int>(n,-1));
  initializeGraph(Graph, n);
  return 0;
}



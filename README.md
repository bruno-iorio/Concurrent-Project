# Concurrent-Project
## Running the code:
There are 3 principal ways to test our code. 
1. main.cpp - To be compiled with ```g++ -std=c++20 -O3 main.cpp -pthread graph.cpp dijkstra.cpp delta_step_static.cpp delta_step_dynamic.cpp -o main``` or by running ```make```. Allows to run an algorithm once, see commandline arguments for parameter specifications. Disclaimer: We used chatGPT writing this part of code, as we deemed it was not core of the project. 
2. benchmarker.py - to be ran with ```python3 benchmarker.py```. Runs main.cpp a number of times with different parameter ranges.
3. alg_compare.cpp - This is the method we recommend for simple testing. Compile with ```g++ -std=c++20 -O3 -pthread graph.cpp dijkstra.cpp delta_step_static.cpp delta_step_dynamic.cpp alg_compare.cpp -o alg_compare``` Disclaimer: We used chatGPT writing this part of code, as we deemed it was not core of the project. The modifiable parameters are as follows, starting at line 60:
  * SRC - source vertex
  * DELTA - delta value
  * THREADS - number of threads
  * name of the roadgraph file you want to supply to the tester
  * afterwards, you can choose if you want to use a roadgraph or a random graph, by uncommenting either line 69 or lines 71 to 76

## Core files:
1. delta_step_static.cpp - implements sequential and parallel static delta stepping with invariant delta. In our testing, only the parallel variant is used, and to simulate the sequential version we set the number of threads to 1.
2. delta_step_dynamic.cpp - corresponding implementations with adaptive delta.
3. dijkstra.cpp - dijkstra for reference.
4. graph.cpp/hpp - function for parsing graphs, and functions for generating random graphs. For one of the functions, wideWeightRandomGraph(), the idea was given by chatGPT (we had a lot of problems debugging our graph generation approach).
## Remarks:
To run the code on roadgraphs, you **need** to create a graphs folder. Within that folder, place .gr files. Compressed .gr files can be found in the road_maps folder on this repo, Otherwise pure .gr files are here: https://drive.google.com/drive/folders/1FwBzXuSA6Vf93w2rZlDwk9j8Adu-oOSI?usp=sharing

## Group Members:
- Bruno Iorio
- Marcel Chwialkowski
- Anca Sfia



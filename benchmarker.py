import subprocess
import re
import pandas as pd
import numpy as np

algo_name_dic = {
    0 : "dijkstra",
    1 : "sds",
    2 : "psds"
}


binary = "./main" 

threads_list = list(range(9,11))
delta_list   = list(range(10000, 20001, 1000))

#options: 0 - dijkstra, 1 - sequential delta stepping, 2 - parallel statting delta stepping        
option = 2
#options: 0 - Random Graph, 1 - RMAT        
option_graph = 0
#how many runs we take per settings -> useful since it varies!
avg = 10

results = []
print("starting tests on " + algo_name_dic[option] + "\n")
for threads in threads_list:
    for delta in delta_list:
        print("starting for", threads, " ", delta, "\n")
        avg_time = 0
        for _ in range(avg):
            cmd = [binary, str(threads), str(delta), str(option)]
            avg_time += float(subprocess.check_output(cmd, text=True))
        avg_time /= avg
        results.append({
            'threads': threads,
            'delta': delta,
            'time_es' : avg_time

        })
        print("done for ", threads, " ", delta, "\n")


#graph_name = "RG" if option_graph == 0 else "RMAT"
#algo_name = algo_name_dic[option]
#results = []
#print("starting tests on " + algo_name_dic[option] + " for " + graph_name + "\n")
#for threads in threads_list:
#    for delta in delta_list:
#        print("starting for", threads, " ", delta, "\n")
#        avg_time = 0
#        for _ in range(avg):
#            cmd = [binary, str(threads), str(delta), str(option), str(option_graph)]
#            avg_time += float(subprocess.check_output(cmd, text=True))
#        avg_time /= avg
#        results.append({
#            'threads': threads,
#            'delta': delta,
#            'time_es' : avg_time
#
#        })
#       print("done for ", threads, " ", delta, "\n")


df = pd.DataFrame(results)
results_name = algo_name + "_" + graph_name + "_" + "results.csv"
df.to_csv('benchmark_results_parallel.csv', index=False)
print("Results saved to " + results_name)

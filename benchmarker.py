import subprocess
import re
import pandas as pd
import numpy as np

algo_name_dic = {
    0 : "dijkstra",
    1 : "dssd",
    2 : "dss",
    3 : "dspd",
    4 : "dsps"
}


binary = "./main" 

threads_list = list(range(1, 10))
delta_list   = list(np.linspace(0.1, 10.0,num= 30))

option = 4
#options: 
#       0 - dijkstra, 
#       1 - delta step sequential dynamic,
#       2 - delta step sequential, 
#       3 - delta step parallel dynamic
#       4 - delta step parallel static

option_graph = 0
#options: 0 - Random Graph, 1 - RMATa

avg = 10
#how many runs we take per settings -> useful since it varies!

graph_name = "RandG" if option_graph == 0 else "RMAT"
algo_name = algo_name_dic[option]
results = []

if option not in [0,1,2,3,4] or option_graph not in [0,1]: 
    exit(1)

print("starting tests on " + algo_name_dic[option] + " for " + graph_name)
if (option in [3,4]):
    for threads in threads_list:
        for delta in delta_list:
            print("starting for", threads, " ", delta)
            cmd = [binary, str(threads), str(delta), str(option), str(option_graph),str(avg)]
            avg_time = float(subprocess.check_output(cmd, text=True))/avg
            results.append({
                'threads': threads,
                'delta': delta,
                'time_es' : avg_time

            })
            print("done for ", threads, " ", delta)
elif option in [1,2]:
    threads = 0
    for delta in delta_list:
        print("starting for ", delta, )
        cmd = [binary, str(threads), str(delta), str(option), str(option_graph),str(avg)]
        avg_time = float(subprocess.check_output(cmd, text=True))/avg
        results.append({
            'delta': delta,
            'time_es' : avg_time
        })
        print("done for ", delta)

elif option == 0:
    threads = 0
    delta = 0
    print("starting for dijkstra")
    cmd = [binary, str(threads), str(delta), str(option), str(option_graph),str(avg)]
    avg_time = float(subprocess.check_output(cmd, text=True))/avg
    
    results.append({
        'time_es' : avg_time
    })
    print("done for dijkstra")

df = pd.DataFrame(results)
results_name = algo_name + "_" + graph_name + "_" + "results.csv"
df.to_csv(results_name, index=False)
print("Results saved to " + results_name)

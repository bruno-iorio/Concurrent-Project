import subprocess
import re
import pandas as pd

binary = "./main" 

threads_list = list(range(8, 9))      
delta_list   = list(range(15000, 20000, 100))   

#options: 0 - dijkstra, 1 - sequential delta stepping, 2 - parallel statting delta stepping        
option = 2   

#how many runs we take per settings -> useful since it varies!
avg = 2

results = []
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

df = pd.DataFrame(results)
df.to_csv('benchmark_results_parallel.csv', index=False)
print("Results saved to benchmark_results_parallel.csv")

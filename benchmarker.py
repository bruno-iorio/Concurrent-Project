import subprocess
import re
import pandas as pd
import numpy as np

test_density = False
default_density = 16
default_delta = 1
default_threads = 10
default_light = 60
default_rebuild = 5

algo_name_dic = {
    0 : "dijkstra",
    1 : "dssd",
    2 : "dss",
    3 : "dspd",
    4 : "dsps"
}


binary = "./main" 
threads_list = list(range(2, 10))
delta_list   = list(np.linspace(0.1, 10.0,num= 30))
density_list = [1 << k for k in range(1,6)]
light_list   = [k for k in range(10,200,20)]
rebuilds_list= [k for k in range(1,10)] 
delta_list3 = [k for k in range(1,10)]

option = 3
#options: 
#       0 - dijkstra, 
#       1 - delta step sequential dynamic,
#       2 - delta step sequential, 
#       3 - delta step parallel dynamic
#       4 - delta step parallel static

option_graph = 0
#options: 0 - Random Graph, 1 - RMATa

avg = 5
#how many runs we take per settings -> useful since it varies!

graph_name = "RandG" if option_graph == 0 else "RMAT"
algo_name = algo_name_dic[option]
results = []

if option not in [0,1,2,3,4] or option_graph not in [0,1]: 
    exit(1)

print("starting tests on " + algo_name_dic[option] + " for " + graph_name)
if (test_density):
    for c in density_list:
        print("starting for desinty", c)
        cmd = [binary, str(default_threads), str(default_delta), str(option), str(option_graph),str(avg),str(c),str(default_light), str(default_rebuild)]
        avg_time = float(subprocess.check_output(cmd, text=True))/avg
        results.append({
                'threads': default_threads,
                'delta': default_delta,
                'time_es' : avg_time

        })
    print("done")
elif (option == 3):
    for threads in threads_list:
        for delta in delta_list3: # because this is integer
            for light in light_list:
                for rebuilds in rebuilds_list:
                    print("starting for", threads, " ", delta)
                    cmd = [binary, str(threads), str(delta), str(option), str(option_graph),str(avg),str(default_density),str(light),str(rebuilds)]
                    avg_time = float(subprocess.check_output(cmd, text=True))/avg

                    results.append({
                        'threads':   threads,
                        'delta':     delta,
                        'rebuilds' : rebuilds,
                        'light' :    light,
                        'time_es' :  avg_time

                    })

                    print("done for ", threads, " ", delta)

elif (option == 4):
    for threads in threads_list:
        for delta in delta_list:
            print("starting for", threads, " ", delta)
            cmd = [binary, str(default_threads), str(default_delta), str(option), str(option_graph),str(avg),str(c),str(default_light), str(default_rebuild)]
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
        cmd = [binary, str(default_threads), str(default_delta), str(option), str(option_graph),str(avg),str(c),str(default_light), str(default_rebuild)]
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
    cmd = [binary, str(threads), str(delta), str(option), str(option_graph),str(avg),str(default_density)]
    avg_time = float(subprocess.check_output(cmd, text=True))/avg
    results.append({
        'time_es' : avg_time
    })
    print("done for dijkstra")



df = pd.DataFrame(results)
results_name = algo_name + "_" + graph_name + "_" + "results.csv"
if (test_density):
    results_name = algo_name + "_" + graph_name + "density_results.csv"
df.to_csv(results_name, index=False)
print("Results saved to " + results_name)

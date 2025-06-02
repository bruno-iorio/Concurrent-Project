import subprocess
import pandas as pd
import numpy as np

# Configuration
binary = "./main"
avg_runs = 2

# Default Parameters
DEFAULTS = {
    'threads': 2,
    'delta': 200,
    'density': 32,
    'light': 60,
    'rebuilds': 5
}

# Options
ALGO_NAMES = {
    0: "dijkstra",
    1: "dss",    # delta step sequential static
    2: "dsps",   # delta step parallel static
    3: "dspd",   # delta step parallel dynamic
}

GRAPH_NAMES = {
    0: "RandG",
    1: "RMAT" # unused unfortunately
}

ROADMAP_NAMES = {
"COL" : "USA-road-d.COL.gr",
"CAL" : "USA-road-d.CAL.gr",
"FLA" : "USA-road-d.FLA.gr",
"LKS" : "USA-road-d.LKS.gr",
"NE"  : "USA-road-d.NE.gr",
"NW"  : "USA-road-d.NW.gr"
}

# Parameter sweeps
threads_list     = list(range(2, 10))
delta_list_long  = list(range(1,300,30))
delta_list_short = list(np.linspace(0.1, 1.0, num=8))
density_list     = [1 << k for k in range(3, 8)] # for density tests
light_list       = list(range(10, 200, 20))
rebuild_list     = list(range(1, 10))

# User options
option_algo = 3        # Set 0-3
option_graph = 0       # Set 0 or 1
option_roadmap = "CAL"

fixed_thread = 2
fixed_delta = 200

test_density = False
test_long_delta = False
test_roadmap = True
test_rebuild_and_light = False
test_fixed_thread = True
test_fixed_delta = True

# Sanity checks
if option_algo not in ALGO_NAMES or option_graph not in GRAPH_NAMES or option_roadmap not in ROADMAP_NAMES:
    raise ValueError("Invalid algorithm or graph option selected")

algo_name = ALGO_NAMES[option_algo]
graph_name = GRAPH_NAMES[option_graph]
graph_name = graph_name if (test_roadmap == False or test_density == True) else option_roadmap
road_file = ROADMAP_NAMES[option_roadmap]
results = []

delta_list = delta_list_long if (test_long_delta) else delta_list_short
delta_list = [DEFAULTS['delta']] if ((not test_roadmap) and test_rebuild_and_light and (option_algo == 3)) else delta_list # this will only test these two parameters

if test_roadmap == False or test_density == True:
    road_file = "0"

if test_fixed_thread == True:
    threads_list = [fixed_thread]

if test_fixed_delta == True:
    delta_list = [fixed_delta]

print(f"Starting benchmark for '{algo_name}' on graph '{graph_name}'\n")
def run_and_get_avg(cmd, avg=avg_runs):
    try:
        output = subprocess.check_output(cmd, text=True)
        return float(output) / avg
    except subprocess.CalledProcessError as e:
        print(f"Error during subprocess: {e}")
        return None

def append_result(**kwargs):
    results.append(kwargs)

# --------------------------- Run Modes ----------------------------

if test_density:
    for density in density_list:
        print(f"Testing density {density}...")
        cmd = [binary, str(DEFAULTS['threads']), str(DEFAULTS['delta']), str(option_algo),
               str(avg_runs), str(density), str(DEFAULTS['light']), str(DEFAULTS['rebuilds']),
               road_file]
        avg_time = run_and_get_avg(cmd)
        if avg_time is not None:
            append_result(threads=DEFAULTS['threads'], delta=DEFAULTS['delta'], density=density, time_es=avg_time)

elif option_algo == 0:  # Dijkstra
    print("Running Dijkstra...")
    cmd = [binary, "0", "0", "0", str(avg_runs), str(DEFAULTS['density']), "0", "0" ,road_file]
    avg_time = run_and_get_avg(cmd)
    if avg_time is not None:
        append_result(time_es=avg_time)

elif option_algo == 1:  # Sequential
    for delta in delta_list:
        print(f"Δ: {delta:.2f}")
        cmd = [binary, "0", str(delta), str(option_algo),
               str(avg_runs), str(DEFAULTS['density']),
               str(0), str(0), road_file]
        avg_time = run_and_get_avg(cmd)
        if avg_time is not None:
            append_result(delta=delta, time_es=avg_time)

elif option_algo == 2:  # Parallel Static
    for threads in threads_list:
        for delta in delta_list:
            print(f"Threads: {threads}, Δ: {delta:.2f}")
            cmd = [binary, str(threads), str(delta), str(option_algo),
                   str(avg_runs), str(DEFAULTS['density']),"0", 
                   "0", road_file]
            avg_time = run_and_get_avg(cmd)
            if avg_time is not None:
                append_result(threads=threads, delta=delta, time_es=avg_time)

elif option_algo == 3:  # Parallel Dynamic
    for threads in threads_list:
        for delta in delta_list:
            if (test_rebuild_and_light):
                for light in light_list:    
                    for rebuild in rebuild_list:
                        print(f"Threads: {threads}, Δ: {delta:.2f}, Light: {light}, Rebuilds: {rebuilds}")
                        cmd = [binary, str(threads), str(delta), str(option_algo),
                               str(avg_runs), str(DEFAULTS['density']),
                               str(light), str(rebuilds), road_file]
                        avg_time = run_and_get_avg(cmd)
                        if avg_time is not None:
                            append_result(threads=threads, delta=delta, light=light,
                                          rebuilds=rebuilds, time_es=avg_time)
            else: 
                print(f"Threads: {threads}, Δ: {delta:.2f}")
                cmd = [binary, str(threads), str(delta), str(option_algo),
                       str(avg_runs), str(DEFAULTS['density']),
                       str(DEFAULTS['light']), str(DEFAULTS['rebuilds']), road_file]
                avg_time = run_and_get_avg(cmd)
                if avg_time is not None:
                    append_result(threads=threads, delta=delta, light=DEFAULTS['light'],
                                  rebuilds=DEFAULTS['rebuilds'], time_es=avg_time)

# -------------------------- Save Results ---------------------------
df = pd.DataFrame(results)
file_suffix = "results"
file_suffix = "rebuild_light_results" if (test_rebuild_and_light and (algo_option==3)) else file_suffix
file_suffix = file_suffix + "_fixedthread" if test_fixed_thread else file_suffix
file_suffix = "density_results" if test_density else file_suffix
csv_name = f"{algo_name}_{graph_name}_{file_suffix}.csv"
df.to_csv(csv_name, index=False)
print(f"\n✅ Results saved to {csv_name}")

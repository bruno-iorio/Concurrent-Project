import pandas as pd
import matplotlib.pyplot as plt

csv_file = 'benchmark_results_parallel.csv'

df = pd.read_csv(csv_file)

plt.figure(figsize=(10, 6))
for threads in sorted(df['threads'].unique()):
    subset = df[df['threads'] == threads]
    plt.plot(subset['delta'], subset['time_es'], label=f"{threads} threads")

plt.xlabel('bucket size (delta)')
plt.ylabel('runtime in microseconds')
plt.title('Parallel delta stepping (static) on the road network in Bay Area')
plt.legend(title='threads')
plt.tight_layout()
plt.show()

"""Generates test data."""

import pandas as pd
import numpy as np

NUM_CELLS = 20
CSV_FILE = 'typical_data.csv'
MARKER_FILE = 'markers.txt'
TSNE_FILE = 'tsne.txt'
CLUSTER_FILE = 'cluster.txt'
cells = pd.Series(['cell_' + str(i) for i in range(1, NUM_CELLS + 1)])
cluster = pd.Series(np.random.randint(1, 4, size=NUM_CELLS))
tSNE_1 = pd.Series(np.random.rand(NUM_CELLS)) + (cluster == 2)
tSNE_2 = pd.Series(np.random.rand(NUM_CELLS)) + (cluster == 3)
gene = pd.DataFrame()
for cluster_index in range(1, 4):
    gene['cluster_' + str(cluster_index)] = (
        # ((cluster == cluster_index) * 10)
        pd.Series(np.random.rand(NUM_CELLS) * 5)
        + ((cluster == cluster_index) * 10)
    )
    gene['not_cluster_' + str(cluster_index)] = (
        ((cluster != cluster_index) * 10)
        # pd.Series(np.random.rand(NUM_CELLS) * 5)
        # + ((cluster != cluster_index) * 10)
    )
    gene['maybe_cluster_' + str(cluster_index)] = (
        pd.Series(np.random.rand(NUM_CELLS) * 10)
        + ((cluster == cluster_index) * 5)
    )

for i in range(1, 9):
    gene['random_' + str(i)] = pd.Series(np.random.rand(NUM_CELLS) * 15)

gene['known_1'] = pd.Series(np.zeros(NUM_CELLS) + 0.0)
gene['known_2'] = pd.Series(np.zeros(NUM_CELLS) + 5.0)
gene['known_3'] = pd.Series(np.zeros(NUM_CELLS) + 0.0)
gene['known_3'][:3] = 2.5
gene['known_3'][3:10] += ((cluster == 1) * 5)

data = pd.DataFrame(
    {'cell': cells, 'cluster': cluster,
     'tSNE_1': tSNE_1, 'tSNE_2':
     tSNE_2}
).join(gene).set_index('cell').rename_axis(None).round(3)
data['cluster'].to_csv(CLUSTER_FILE)
print("Saved to " + CLUSTER_FILE)
data[['tSNE_1', 'tSNE_2']].to_csv(TSNE_FILE, header=None)
print("Saved to " + TSNE_FILE)
data[data.columns[3:]].to_csv(MARKER_FILE)
print("Saved to " + MARKER_FILE)
for gene in data.columns[3:]:
    data[gene + "_c"] = -data[gene]

data.to_csv(CSV_FILE)
print("Saved to " + CSV_FILE)

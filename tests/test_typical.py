"""A basic pytest test case for hgmd.


Generates a "typical"-looking dataset with following specifications:
    - 20 cells
    - 20 genes
    - 3 clusters
    - 1 gene clearly (but not absolutely) marks each cluster
    - 1 gene's complement does the same
    - 1 gene roughly marks each cluster
    - 1 gene has no expression
    - 1 gene has near uniform expression
    - 5 genes with known values
    - 6 genes with random values
Verifies that all functions in hgmd.py return what they should return. Does not
check for exceptions (yet!). Does not check for PDF correctness (yet!!).
"""
# TODO: test fringe cases: weird data, exception-raising (update docstring!).
# TODO: PDFs (also update docstring!).

import pandas as pd
import numpy as np

from hgmd import hgmd as hg


def setup_module():
    """Generates CSVs of all test data."""

    # TEST_DATA_FILE contains cluster, tSNE 1/2, and known gene data
    TEST_DATA_FILE = 'typical_data.csv'
    MARKER_FILE = 'markers.txt'
    TSNE_FILE = 'tsne.txt'
    CLUSTER_FILE = 'cluster.txt'
    NUM_CELLS = 20

    data = pd.read_csv(TEST_DATA_FILE)
    cluster = data['cluster']
    gene = pd.DataFrame()

    for cluster_index in range(1, 4):
        gene['cluster_' + str(cluster_index)] = (
            pd.Series(np.random.rand(NUM_CELLS))
            + ((cluster == cluster_index) * 10)
        )
        gene['not_cluster_' + str(cluster_index)] = (
            pd.Series(np.random.rand(NUM_CELLS))
            + ((cluster != cluster_index) * 10)
        )
        gene['maybe_cluster_' + str(cluster_index)] = (
            pd.Series(np.random.rand(NUM_CELLS) * 5)
            + ((cluster == cluster_index) * 5)
        )

    for i in range(1, 7):
        gene['random_' + str(i)] = pd.Series(
            np.random.rand(NUM_CELLS) * 15
        )

    data = data.join(gene).set_index('cell').rename_axis(None)
    data['cluster'].to_csv(CLUSTER_FILE)
    data[['tSNE_1', 'tSNE_2']].to_csv(TSNE_FILE)
    data[data.columns[3:]].to_csv(MARKER_FILE)


class TestTypical:
    def test_func(self):
        x = 3
        assert x is 4

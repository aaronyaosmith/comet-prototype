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
import pytest

from hgmd import hgmd as hg


# TEST_DATA_FILE contains cluster, tSNE 1/2, and known gene data
FOLDER = 'typical_data/'
TEST_DATA_FILE = FOLDER + 'typical_data.csv'
MARKER_FILE = FOLDER + 'markers.txt'
TSNE_FILE = FOLDER + 'tsne.txt'
CLUSTER_FILE = FOLDER + 'cluster.txt'
SINGLETON_OUTPUT = FOLDER + 'singleton/'
PAIR_OUTPUT = FOLDER + 'pair/'
TP_TN_OUTPUT = FOLDER + 'TP_TN/'
NUM_CELLS = 20


@pytest.fixture(scope='module')
def csv_read_data():
    return pd.read_csv(TEST_DATA_FILE, index_col=0).rename_axis(None)


def assert_frame_equal(df1, df2):
    """Equality for frames containing floats"""
    df1.fillna(0)
    df2.fillna(0)
    return pd.testing.assert_frame_equal(df1, df2, check_dtype=False)


class TestTypical:
    def test_get_cell_data(self, csv_read_data):
        func_data = hg.get_cell_data(
            marker_path=MARKER_FILE,
            tsne_path=TSNE_FILE,
            cluster_path=CLUSTER_FILE
        )
        assert_frame_equal(csv_read_data, func_data)

    def test_singleton_test(self, csv_read_data):
        for cluster in csv_read_data['cluster'].unique():
            func_data = hg.singleton_test(csv_read_data, cluster, 0, NUM_CELLS)
            path = SINGLETON_OUTPUT + "cluster_" + str(cluster) + ".csv"
            read_data = pd.read_csv(path, index_col=0).rename_axis(None)
            assert_frame_equal(read_data, func_data)

    def test_pair_test(self, csv_read_data):
        for cluster in csv_read_data['cluster'].unique():
            singleton_path = (
                SINGLETON_OUTPUT + "cluster_" + str(cluster) + ".csv"
            )
            singleton = pd.read_csv(
                singleton_path, index_col=0
            ).rename_axis(None)
            func_data = hg.pair_test(csv_read_data, singleton, cluster)
            func_data.to_csv(PAIR_OUTPUT + "cluster_" + str(cluster) + ".csv")

    def test_find_TP_TN(self, csv_read_data):
        for cluster in csv_read_data['cluster'].unique():
            singleton_path = (
                SINGLETON_OUTPUT + "cluster_" + str(cluster) + ".csv"
            )
            singleton = pd.read_csv(
                singleton_path, index_col=0
            ).rename_axis(None)
            pair_path = (
                PAIR_OUTPUT + "cluster_" + str(cluster) + ".csv"
            )
            pair = pd.read_csv(
                pair_path, index_col=0
            ).rename_axis(None)
            singleton, pair = hg.find_TP_TN(
                csv_read_data, singleton, pair, cluster
            )
            singleton.to_csv(
                TP_TN_OUTPUT + "singleton_cluster_" + str(cluster) + ".csv"
            )
            pair.to_csv(TP_TN_OUTPUT + "pair_cluster_" + str(cluster) + ".csv")

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

The OUTPUT_FOLDER stores "correct" output to be checked
against. TEST_OUTPUT_FOLDER stores the output of the test runs.
"""
# TODO: test fringe cases: weird data, exception-raising (update docstring!).
# TODO: PDFs (also update docstring!).

import pandas as pd
import pytest

from hgmd import hgmd as md


# TEST_DATA_FILE contains cluster, tSNE 1/2, and known gene data
INPUT_FOLDER = '../../data/testing/typical_data/input/'
OUTPUT_FOLDER = '../../data/testing/typical_data/output/'
TEST_OUTPUT_FOLDER = '../../data/testing/typical_data/test_output/'
TEST_DATA_FILE = '../../data/testing/typical_data/typical_data.csv'
MARKER_FILE = INPUT_FOLDER + 'markers.txt'
TSNE_FILE = INPUT_FOLDER + 'tsne.txt'
CLUSTER_FILE = INPUT_FOLDER + 'cluster.txt'


@pytest.fixture(scope='module')
def csv_read_data():
    return pd.read_csv(TEST_DATA_FILE, index_col=0).rename_axis(None)


def assert_frame_equal(df1, df2):
    """Equality for frames containing floats"""
    df1.fillna(0)
    df2.fillna(0)
    return pd.testing.assert_frame_equal(
        df1, df2, check_dtype=False
    )


class TestTypical:
    @staticmethod
    def singleton_data(cluster):
        path = (
            OUTPUT_FOLDER + 'cluster_' + str(cluster)
            + '/singleton_data.csv'
        )
        expected = pd.read_csv(path, index_col=0).rename_axis(None)
        expected = expected[[
            'gene', 'HG_stat', 'mHG_pval', 'mHG_cutoff_index',
            'mHG_cutoff_value',
            't_stat', 't_pval', 'HG_rank', 'rank', 't_rank',
            'sequential_rank'
        ]]
        return expected

    @staticmethod
    def pair_data(cluster):
        path = (
            OUTPUT_FOLDER + 'cluster_' + str(cluster)
            + '/pair_data.csv'
        )
        expected = pd.read_csv(path, index_col=0).rename_axis(None)
        expected = expected[[
            'gene', 'gene_B', 'HG_stat', 'mHG_pval', 'mHG_cutoff_index',
            'mHG_cutoff_value',
            't_stat', 't_pval', 'HG_rank', 'rank', 't_rank',
            'sequential_rank'
        ]]
        return expected

    @staticmethod
    def TP_TN_singleton_data(cluster):
        path = (
            OUTPUT_FOLDER + 'cluster_' + str(cluster)
            + '/singleton_data.csv'
        )
        expected = pd.read_csv(path, index_col=0).rename_axis(None)
        expected = expected[[
            'gene', 'HG_stat', 'mHG_pval', 'mHG_cutoff_index',
            'mHG_cutoff_value',
            't_stat', 't_pval', 'HG_rank', 'rank', 't_rank',
            'sequential_rank',
            'true_negative', 'true_positive'
        ]]
        return expected

    @staticmethod
    def TP_TN_data(cluster):
        path = (
            OUTPUT_FOLDER + 'cluster_' + str(cluster)
            + '/pair_data.csv'
        )
        expected = pd.read_csv(path, index_col=0).rename_axis(None)
        expected = expected[[
            'gene', 'gene_B', 'HG_stat', 'mHG_pval', 'mHG_cutoff_index',
            'mHG_cutoff_value',
            't_stat', 't_pval', 'HG_rank', 'rank', 't_rank',
            'sequential_rank',
            'true_negative', 'true_positive'
        ]]
        return expected

    @staticmethod
    def TP_TN_weighted_singleton_data(cluster):
        path = (
            OUTPUT_FOLDER + 'cluster_' + str(cluster)
            + '/singleton_data.csv'
        )
        expected = pd.read_csv(path, index_col=0).rename_axis(None)
        return expected

    @staticmethod
    def TP_TN_weighted_data(cluster):
        path = (
            OUTPUT_FOLDER + 'cluster_' + str(cluster)
            + '/pair_data.csv'
        )
        expected = pd.read_csv(path, index_col=0).rename_axis(None)
        return expected

    def test_get_cell_data(self, csv_read_data):
        func_data = md.get_cell_data(
            marker_path=MARKER_FILE,
            tsne_path=TSNE_FILE,
            cluster_path=CLUSTER_FILE
        )
        assert_frame_equal(csv_read_data, func_data)

    def test_singleton_test(self, csv_read_data):
        for cluster in csv_read_data['cluster'].unique():
            func_data = md.singleton_test(csv_read_data, cluster)
            func_data = func_data.reset_index(drop=True)
            expected = self.singleton_data(cluster)
            func_data.to_csv(
                TEST_OUTPUT_FOLDER + 'cluster_' + str(cluster)
                + '_singleton.csv'
            )
            assert_frame_equal(expected, func_data)

    def test_pair_test(self, csv_read_data):
        for cluster in csv_read_data['cluster'].unique():
            singleton = self.singleton_data(cluster)
            func_data = md.pair_test(csv_read_data, singleton, cluster)
            func_data = func_data.reset_index(drop=True)
            expected = self.pair_data(cluster)
            func_data.to_csv(
                TEST_OUTPUT_FOLDER + 'cluster_' + str(cluster) + '_pair.csv'
            )
            assert_frame_equal(expected, func_data)

    def test_find_TP_TN(self, csv_read_data):
        for cluster in csv_read_data['cluster'].unique():
            singleton = self.singleton_data(cluster)
            pair = self.pair_data(cluster)
            singleton, pair = md.find_TP_TN(
                csv_read_data, singleton, pair, cluster
            )
            singleton_expected = self.TP_TN_singleton_data(cluster)
            pair_expected = self.TP_TN_data(cluster)
            singleton.to_csv(
                TEST_OUTPUT_FOLDER + 'cluster_' + str(cluster)
                + '_TP_TN_singleton.csv'
            )
            pair.to_csv(
                TEST_OUTPUT_FOLDER + 'cluster_' + str(cluster)
                + '_TP_TN_pair.csv'
            )
            assert_frame_equal(singleton_expected, singleton)
            assert_frame_equal(pair_expected, pair)

    def test_find_weighted_TP_TN(self, csv_read_data):
        for cluster in csv_read_data['cluster'].unique():
            singleton = self.TP_TN_weighted_singleton_data(cluster)
            pair = self.TP_TN_weighted_data(cluster)
            singleton, pair = md.find_TP_TN(
                csv_read_data, singleton, pair, cluster
            )
            singleton_expected = self.TP_TN_weighted_singleton_data(cluster)
            pair_expected = self.TP_TN_weighted_data(cluster)
            singleton.to_csv(
                TEST_OUTPUT_FOLDER + 'cluster_' + str(cluster)
                + '_TP_TN_weighted_singleton.csv'
            )
            pair.to_csv(
                TEST_OUTPUT_FOLDER + 'cluster_' + str(cluster)
                + '_TP_TN_weighted_pair.csv'
            )
            assert_frame_equal(singleton_expected, singleton)
            assert_frame_equal(pair_expected, pair)

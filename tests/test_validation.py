"""Tests input validation with various bad input."""

import pytest

from hgmd import hgmd as md


INPUT_FOLDER = '../../data/testing/validation_data/'


def test_missing_file():
    with pytest.raises(OSError) as excinfo:
        marker_path = INPUT_FOLDER + 'WRONG_MARKERS.txt'
        tsne_path = INPUT_FOLDER + 'tsne.txt'
        cluster_path = INPUT_FOLDER + 'cluster.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'WRONG_MARKERS.txt' in str(excinfo.value)
    with pytest.raises(OSError) as excinfo:
        marker_path = INPUT_FOLDER + 'markers.txt'
        tsne_path = INPUT_FOLDER + 'WRONG_TSNE.txt'
        cluster_path = INPUT_FOLDER + 'cluster.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'WRONG_TSNE.txt' in str(excinfo.value)
    with pytest.raises(OSError) as excinfo:
        marker_path = INPUT_FOLDER + 'markers.txt'
        tsne_path = INPUT_FOLDER + 'tsne.txt'
        cluster_path = INPUT_FOLDER + 'WRONG_CLUSTER.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'WRONG_CLUSTER.txt' in str(excinfo.value)


def test_invalid_path():
    with pytest.raises(OSError) as excinfo:
        marker_path = 'hello'
        tsne_path = INPUT_FOLDER + 'tsne.txt'
        cluster_path = INPUT_FOLDER + 'cluster.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'hello' in str(excinfo.value)
    with pytest.raises(OSError) as excinfo:
        marker_path = INPUT_FOLDER + 'markers.txt'
        tsne_path = INPUT_FOLDER + 'DUMMY/'
        cluster_path = INPUT_FOLDER + 'cluster.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'DUMMY/' in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        marker_path = INPUT_FOLDER + 'markers.txt'
        tsne_path = INPUT_FOLDER + 'tsne.txt'
        cluster_path = 0
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'int' in str(excinfo.value)


def test_invalid_format():
    with pytest.raises(ValueError) as excinfo:
        marker_path = INPUT_FOLDER + 'markers_spaces.csv'
        tsne_path = INPUT_FOLDER + 'tsne.txt'
        cluster_path = INPUT_FOLDER + 'cluster.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'comma' in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        marker_path = INPUT_FOLDER + 'markers.txt'
        tsne_path = INPUT_FOLDER + 'tsne_wrong_column_order.csv'
        cluster_path = INPUT_FOLDER + 'cluster.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    with pytest.raises(ValueError) as excinfo:
        marker_path = INPUT_FOLDER + 'markers.txt'
        tsne_path = INPUT_FOLDER + 'tsne.txt'
        cluster_path = INPUT_FOLDER + 'cluster_extra_column.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'has too many columns' in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        marker_path = INPUT_FOLDER + 'markers.txt'
        tsne_path = INPUT_FOLDER + 'tsne.txt'
        cluster_path = INPUT_FOLDER + 'cluster_missing_column.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'missing column' in str(excinfo.value)
    # TODO: validate nonsense: corrupt file, messed up header, etc.
    """
    with pytest.raises(ValueError) as excinfo:
        marker_path = INPUT_FOLDER + 'NONSENSE.txt'
        tsne_path = INPUT_FOLDER + 'tsne.txt'
        cluster_path = INPUT_FOLDER + 'cluster_missing_column.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    """


def test_invalid_data():
    with pytest.raises(ValueError) as excinfo:
        marker_path = INPUT_FOLDER + 'markers_missing_data.csv'
        tsne_path = INPUT_FOLDER + 'tsne.txt'
        cluster_path = INPUT_FOLDER + 'cluster.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'missing' in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        marker_path = INPUT_FOLDER + 'markers.txt'
        tsne_path = INPUT_FOLDER + 'tsne_cell_mismatch.csv'
        cluster_path = INPUT_FOLDER + 'cluster.txt'
        md.get_cell_data(
            marker_path, tsne_path, cluster_path
        )
    assert 'inconsistent' in str(excinfo.value)


def test_warnings():
    # TODO: warn if weird but valid input
    pass

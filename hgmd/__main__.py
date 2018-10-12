import os
import argparse
import datetime

import pandas as pd

from . import hgmd
from docs.source import conf


def main():
    """Hypergeometric marker detection. Finds markers identifying a cluster.

    Reads in data from single-cell RNA sequencing. Data is in the form of 3
    CSVs: gene expression data by gene by cell, 2-D tSNE data by cell, and the
    clusters of interest by cell. Creates a list of genes and a list of gene
    pairs (including complements), ranked by hypergeometric and t-test
    significance. The highest ranked marker genes generally best identify the
    cluster of interest. Saves these lists to CSV and creates gene expression
    visualizations.
    """
    # TODO: more precise description

    start_dt = datetime.datetime.now()
    print("Started on " + str(start_dt.isoformat()))

    parser = argparse.ArgumentParser(
        description=("Hypergeometric marker detection. Finds markers "
                     "identifying a cluster.")
    )
    parser.add_argument(
        'input_path', type=str,
        help=("the input directory containing markers.txt, tsne.txt, and "
              "cluster.txt.")
    )
    parser.add_argument(
        'output_path', type=str,
        help="the output directory where output files should go"
    )
    parser.add_argument(
        '-X', nargs='?', default=None,
        help="X argument for XL-mHG"
    )
    parser.add_argument(
        '-L', nargs='?', default=None,
        help="L argument for XL-mHG"
    )
    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    X = args.X
    L = args.L
    marker_file = 'markers.txt'
    tsne_file = 'tsne.txt'
    cluster_file = 'cluster.txt'

    # TODO: gene pairs with expression ratio within the cluster of interest
    # under [min_exp_ratio] were ignored in hypergeometric testing. This
    # functionality is currently unimplemented.

    # min_exp_ratio = 0.4

    csv_path = output_path + 'data/'
    vis_path = output_path + 'vis/'
    pickle_path = output_path + '_pickles/'
    os.makedirs(csv_path, exist_ok=True)
    os.makedirs(vis_path, exist_ok=True)
    os.makedirs(pickle_path, exist_ok=True)

    if X is not None:
        X = int(X)
        print("Set X to " + str(X) + ".")
    if L is not None:
        L = int(L)
        print("Set L to " + str(L) + ".")
    print("Reading data...")
    cls_ser = pd.read_csv(
        input_path + cluster_file,
        index_col=0, names=['cell', 'cluster'], squeeze=True
    )

    tsne = pd.read_csv(
        input_path + tsne_file,
        index_col=0, names=['cell', 'tSNE_1', 'tSNE_2']
    )
    no_complement_marker_exp = pd.read_csv(
        input_path + marker_file, index_col=0
    ).rename_axis('cell')

    print("Generating complement data...")
    marker_exp = hgmd.add_complements(no_complement_marker_exp)

    # Enumerate clusters and process each individually in its own folder.
    clusters = cls_ser.unique()
    clusters.sort()
    for cls in clusters:
        print('Processing cluster ' + str(cls) + '...')
        print('Running t test on singletons...')
        t_test = hgmd.batch_t(marker_exp, cls_ser, cls)
        print('Running XL-mHG on singletons...')
        xlmhg = hgmd.batch_xlmhg(marker_exp, cls_ser, cls, X=X, L=L)
        # We need to slide the cutoff indices before using them,
        # to be sure they can be used in the real world. See hgmd.mhg_slide()
        cutoff_value = hgmd.mhg_cutoff_value(
            marker_exp, xlmhg[['gene', 'mHG_cutoff']]
        )
        xlmhg = xlmhg[['gene', 'HG_stat', 'mHG_pval']].merge(
            hgmd.mhg_slide(marker_exp, cutoff_value), on='gene'
        )
        # Update cutoff_value after sliding
        cutoff_value = pd.Series(
            xlmhg['cutoff_val'].values, index=xlmhg['gene'])
        print('Creating discrete expression matrix...')
        discrete_exp = hgmd.discrete_exp(marker_exp, cutoff_value)
        print('Finding pair expression matrix...')
        (
            gene_map, in_cls_count, pop_count,
            in_cls_product, total_product, upper_tri_indices
        ) = hgmd.pair_product(discrete_exp, cls_ser, cls)
        print('Running hypergeometric test on pairs...')
        pair = hgmd.pair_hg(
            gene_map, in_cls_count, pop_count,
            in_cls_product, total_product, upper_tri_indices
        )
        print('Finding simple true positives/negatives for singletons...')
        sing_tp_tn = hgmd.tp_tn(discrete_exp, cls_ser, cls)
        print('Finding simple true positives/negatives for pairs...')
        pair_tp_tn = hgmd.pair_tp_tn(
            gene_map, in_cls_count, pop_count,
            in_cls_product, total_product, upper_tri_indices
        )

        # Save TP/TN values to be used for non-cluster-specific things
        print('Pickling data for later...')
        sing_tp_tn.to_pickle(pickle_path + 'sing_tp_tn_' + str(cls))
        pair_tp_tn.to_pickle(pickle_path + 'pair_tp_tn_' + str(cls))

        # Export to CSV for user
        print('Exporting cluster ' + str(cls) + ' output to CSV...')
        sing_output = xlmhg\
            .merge(t_test, on='gene')\
            .merge(sing_tp_tn, on='gene')\
            .set_index('gene')\
            .sort_values(by='HG_stat', ascending=True)
        sing_output['rank'] = sing_output.reset_index().index + 1
        sing_output.to_csv(
            csv_path + '/cluster_' + str(cls) + '_singleton.csv'
        )
        sing_stripped = sing_output[
            ['HG_stat', 'TP', 'TN']
        ].reset_index().rename(index=str, columns={'gene': 'gene_1'})
        pair_output = pair\
            .merge(pair_tp_tn, on=['gene_1', 'gene_2'], how='left')\
            .append(sing_stripped, sort=False, ignore_index=True)\
            .sort_values(by='HG_stat', ascending=True)
        pair_output['rank'] = pair_output.reset_index().index + 1
        pair_output.to_csv(
            csv_path + '/cluster_' + str(cls) + '_pair.csv'
        )

    # Add text file to keep track of everything
    end_dt = datetime.datetime.now()
    print("Ended on " + end_dt.isoformat())
    metadata = open(output_path + 'metadata.txt', 'w')
    metadata.write("Started: " + start_dt.isoformat())
    metadata.write("\nEnded: " + end_dt.isoformat())
    metadata.write("\nElapsed: " + str(end_dt - start_dt))
    metadata.write("\nGenerated by COMET version " + conf.version)


if __name__ == '__main__':
    main()

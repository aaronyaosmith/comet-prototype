import os

from . import hgmd as md


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

    # TODO: get command line input instead of assigning directly
    input_path = '/home/aaron/Documents/hgmd/data/input/'
    output_path = '/home/aaron/Documents/hgmd/data/output/'
    marker_file = '20180327_col_batch_1_2_3_tpm_combined_marker_mat.txt'
    tsne_file = '20180327_col_batch_1_2_3_CCA_corrected_1_10_tSNE.txt'
    cluster_file = '20180327_col_batch_1_2_3_CCA_corrected_1_10_clusters.txt'
    X = 3
    L = 750
    min_exp_ratio = 0.4
    plot_pages = 60
    plot_genes = 20

    print("Reading data...")
    cell_data = md.get_cell_data(
        marker_path=(input_path + marker_file),
        tsne_path=(input_path + tsne_file),
        cluster_path=(input_path + cluster_file)
    )

    # Enumerate clusters and process each individually in its own folder.
    # pair_data also contains singleton data, but singleton is just
    # singleton.
    clusters = cell_data['cluster'].unique()
    clusters.sort()
    for cluster in clusters:
        print("Processing cluster " + str(cluster) + "...")
        cluster_path = output_path + "_cluster_" + str(cluster)
        os.makedirs(cluster_path, exist_ok=True)
        print("Testing singletons...")
        singleton_data = md.singleton_test(cell_data, cluster, X, L)
        print("Testing pairs...")
        pair_data = md.pair_test(
            cell_data, singleton_data, cluster, min_exp_ratio
        )
        print("Calculating true positive/negative rates...")
        md.find_TP_TN(cell_data, singleton_data, pair_data, cluster)
        print("Saving to CSV...")
        singleton_data.to_csv(cluster_path + "singleton_data.csv")
        pair_data.to_csv(cluster_path + "pair_data.csv")
        print("Done.")
        print("Plotting true positive/negative rates...")
        md.make_TP_TN_plots(
            cell_data, singleton_data, pair_data, plot_genes,
            pair_path=(cluster_path + "TP_TN_plot.pdf"),
            singleton_path=(cluster_path + "singleton_TP_TN_plot.pdf")
        )
        print("Done.")
        print("Plotting discrete expression...")
        md.make_discrete_plots(
            cell_data, singleton_data, pair_data, plot_pages,
            path=(cluster_path + "discrete_plots.pdf"),
        )
        print("Done.")
        print("Plotting continuous expresssion...")
        md.make_combined_plots(
            cell_data, singleton_data, pair_data, plot_genes,
            pair_path=(cluster_path + "combined_plot.pdf"),
            singleton_path=(cluster_path + "singleton_combined_plot.pdf")
        )
        print("Done.")

    print("All set!!! Enjoy your PDFs and CSVs.")


if __name__ == '__main__':
    main()

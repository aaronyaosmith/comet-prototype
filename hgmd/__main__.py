import os
import argparse

from . import hgmd as md
import sys
import multiprocessing 
import time
import math


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
    # TODO: get X and L as CLI args

    # TODO: get command line input instead of assigning directly
    parser = argparse.ArgumentParser(
        description=("Hypergeometric marker detection. Finds markers "
                     "identifying a cluster.")
    )
    parser.add_argument(
        'marker', type=str,
        help=("Marker file input")
    )
    parser.add_argument(
        'tsne', type=str,
        help=("tsne file input")
    )
    parser.add_argument(
        'cluster', type=str,
        help=("Cluster file input")
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
    output_path = args.output_path
    X = args.X
    L = args.L
    marker_file = args.marker
    tsne_file = args.tsne
    cluster_file = args.cluster
    min_exp_ratio = 0.4
    plot_pages = 15
    plot_genes = 15

    if X is not None:
        X = int(X)
        print("Set X to " + str(X) + ".")
    if L is not None:
        L = int(L)
        print("Set L to " + str(L) + ".")
    print("Reading data...")
    cell_data = md.get_cell_data(
        marker_path=(marker_file),
        tsne_path=(tsne_file),
        cluster_path=(cluster_file)
    )

    # Enumerate clusters and process each individually in its own folder.
    # pair_data also contains singleton data, but singleton is just
    # singleton.
    clusters = cell_data['cluster'].unique()
    clusters.sort()
    print('BEEP BOOP')
    print(clusters)
    start_time = time.time()
    
    #cores is number of simultaneous threads you want to run, can be set at will
    cores = 1
    # if core number is bigger than number of clusters, just set it equal to number of clusters
    if cores > len(clusters):
        cores = len(clusters)
    #below loops allow for splitting the job based on core choice
    group_num  = math.ceil((len(clusters) / cores ))
    for element in range(group_num):
        new_clusters = clusters[:cores]
        print(new_clusters)
        jobs = []
        #this loop spawns the workers and runs the code for each assigned.
        #workers assigned based on the new_clusters list which is the old clusters
        #split up based on core number e.g.
        #clusters = [1 2 3 4 5 6] & cores = 4 --> new_clusters = [1 2 3 4], new_clusters = [5 6]
        for cluster in new_clusters:
            p = multiprocessing.Process(target=md.process,
                args=(cluster,cell_data,X,L,min_exp_ratio,plot_pages,plot_genes,output_path))
            jobs.append(p)
            p.start()
        p.join()
        
        new_clusters = []
        clusters = clusters[cores:len(clusters)]
        print(clusters)

    end_time = time.time()

    print('Took ' + str(end_time-start_time) + ' seconds')
    print('Which is ' + str( (end_time-start_time)/60 ) + ' minutes')

    print("All set!!! Enjoy your PDFs and CSVs.")


if __name__ == '__main__':
    main()

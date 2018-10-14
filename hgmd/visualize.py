import pandas as pd
import numpy as np
import xlmhg as hg
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages

"""
Set of modularized visualization functions for COMET; producing graphs in PDFs.
"""

"""
Goals: Separate singleton and pair graphs for:

    - TP/TN with gene names

    - Side-by-side comparison between discrete and continuous expression

For singletons and pairs combined:

    - Side-by-side comparison of single gene discrete expression and resulting
      pair discrete expression.
"""


CMAP_CONTINUOUS = cm.get_cmap('nipy_spectral')
CMAP_DISCRETE = cm.get_cmap('bwr')


def make_titles(genes):
    """Vectorized version of make_title(). Takes DataFrame as input.

    genes should have columns 'gene_1', 'gene_2', 'rank', 'cutoff_val'.
    """

    vmt = np.vectorize(make_title)
    return vmt(
        genes['gene_1'].values, genes['rank'].values,
        genes['gene_2'].values, genes['cutoff_val'].values
    )


def make_title(gene_1, rank, gene_2=None, cutoff_val=None):
    """Makes a plot title for a gene or gene pair.

    Formatting: for pairs, 'rank $rank: $gene_1+$gene_2', and for singletons,
    'rank $rank: $gene_1 $cutoff_val'.  gene_2 should be None for singletons,
    and cutoff_val should be None for pairs.

    :param genes: A DataFrame with columns 'gene_1', 'gene_2', 'rank',
        'cutoff_val'.

    :returns: A list of strings containing the titles, with indices
              corresponding to the input DataFrame.

    :rtype: string list
    """

    if gene_2 is not None and cutoff_val is not None:
        raise Exception("Invalid argument")

    if gene_2 is None:
        return (
            "rank " + str(rank) + ": " + str(gene_1) + " " + str(cutoff_val)
        )
    elif cutoff_val is None:
        return ("rank " + str(rank) + ": " + str(gene_1) + "+" + str(gene_2))
    else:
        raise Exception("Invalid argument")


def make_discrete_plots(tsne, discrete_exp, plot_genes, path):
    """Plots discrete gene expression of paired genes to PDF.

    For each gene pair listed in plot_genes, make three scatterplots.  First, a
    plot showing joint expression.  Then, two plots showing singleton
    expression for each genes.  If a single gene is passed, plot only its
    expression, and make two blank plots.  Save each gene/gene pair as a PDF
    page, then save to path.

    :param tsne: A DataFrame with 'cell', 'tSNE_1', and 'tSNE_2' columns.
    :param discrete_exp: A DataFrame whose rows are cell identifiers, columns
        are gene identifiers, and values are boolean values representing gene
        expression.
    :param plot_genes: A list of tuples, each element containing a tuple
        containing titles of the three graphs, the first gene name, and the
        second gene name. For singletons, the 3rd item is None.
    :param path: The path to which the PDF will be saved.

    :returns: Nothing.
    """

    def make_plot(ax, title, coords, cmap):
        """Make a single graph on ax with given specs."""
        ax.set_title(title)
        ax.set_xlabel('tSNE_1')
        ax.set_ylabel('tSNE_2')
        sc = ax.scatter(
            x=coords[0],
            y=coords[1],
            c=coords[2],
            s=2
        )
        plt.colorbar(sc, ax=ax)

    def make_pair_discrete_page(fig, ax_triple, titles, gene_1, gene_2):
        """Make page with three discrete plots given titles and genes."""
        coords_df = tsne.merge(discrete_exp[[gene_1, gene_2]], on='cell')
        coords_df['pair'] = coords_df[gene_1] * coords_df[gene_2]

        for (graph_index, z_label) in ((0, 'pair'), (1, gene_1), (2, gene_2)):
            make_plot(
                ax=ax_triple[graph_index],
                title=titles[graph_index],
                coords=(
                    coords_df['tSNE_1'].values,
                    coords_df['tSNE_2'].values,
                    coords_df[z_label].values
                ),
                cmap=CMAP_DISCRETE
            )

    def make_single_discrete_page(fig, ax_triple, title, gene):
        """Make page with one discrete plot given title and gene"""
        coords_df = tsne.merge(discrete_exp[gene], on='cell')
        make_plot(
            ax=ax_triple[0],
            title=title[0],
            coords=(
                coords_df['tSNE_1'].values,
                coords_df['tSNE_2'].values,
                coords_df[gene].values
            ),
            cmap=CMAP_DISCRETE
        )

    pdf = PdfPages(path)
    for plot_gene in plot_genes:
        fig, ax_triple = plt.subplots(ncols=3, figsize=(15, 5))
        if len(plot_gene) == 3:
            make_pair_discrete_page(
                fig=fig, ax_triple=ax_triple,
                titles=plot_gene[0],
                gene_1=plot_gene[1],
                gene_2=plot_gene[2]
            )
        elif len(plot_gene) == 2:
            make_single_discrete_page(
                fig=fig, ax_triple=ax_triple,
                title=plot_gene[0][0],
                gene=plot_gene[1]
            )
        else:
            raise Exception("Invalid argument")
        pdf.savefig(fig)
        plt.close(fig)


def make_combined_plots(tsne, discrete_exp, marker_exp, plot_genes, path):
    """Plots discrete alongside continuous expression to PDF.

    For each gene/gene pair listed in plot_genes, make two scatterplots: a plot
    showing discrete expression, and a plot showing expression on a color
    spectrum.  For gene pairs, make these two plots separately for each gene.
    Save each gene/gene pair as a PDF page, then save to path.

    :param tsne: A DataFrame with 'cell', 'tSNE_1', and 'tSNE_2' columns.
    :param discrete_exp: A DataFrame whose rows are cell identifiers, columns
        are gene identifiers, and values are boolean values representing gene
        expression.
    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param plot_genes: An array whose elements are either single gene names, or
        tuples containing two gene names.
    :param path: The path to which the PDF will be saved.

    :returns: Nothing.
    """
    raise Exception("Unimplemented")


def make_TP_TN_plots(plot_genes, sing_tp_tn, pair_tp_tn, path):
    """Plots TP/TN rates of genes/pairs to PDF.

    For each gene/gene pair listed in plot_genes, plot their TP/TN rate on a
    scatterplot, labeling the point with the gene/gene pair name.  When done,
    output this scatterplot to PDF and save to path.

    :param plot_genes: An array whose elements are either single gene names, or
        tuples containing two gene names.
    :param sing_tp_tn: A DataFrame with 'gene', 'TP', and 'TN' columns.
    :param pair_tp_tn: A DataFrame with 'gene_1', 'gene_2', 'TP', and 'TN'
        columns.
    :param path: The path to which the PDF will be saved.

    :returns: Nothing.
    """
    raise Exception("Unimplemented")


"""
def make_discrete_plots(cells, singleton, pair, plot_pages, path):
    Creates plots of discrete expression and saves them to pdf.
    Uses output of singleton_test and pair_test to create plots of discrete
    gene expression. For gene pairs, plot each individual gene as well,
    comparing their plots to the pair plot. Plots most significant genes/pairs
    first.
    Args:
        cells: A DataFrame with format matching those returned by
            get_cell_data. Row values are cell identifiers, columns are first
            cluster identifier, then tSNE_1 and tSNE_2, then gene names.
        singleton: A DataFrame with format matching those returned by
            singleton_test.
        pair: A DataFrame with format matching those returned by pair_test.
        plot_pages: The maximum number of pages to plot. Each gene or gene pair
            corresponds to one page.
        path: Save the plot pdf here.
    Returns:
        Nothing.
    Raises:
        ValueError: cells, singleton, or pair is in an incorrect format,
            plot_pages is less than 1.

    exp = get_discrete_exp(cells, singleton)
    with PdfPages(path) as pdf:
        for i in range(0, plot_pages):
            print("Plotting discrete plot "+str(i+1)+" of "+str(plot_pages))
            gene_A = pair['gene'].iloc[i]
            gene_B = pair['gene_B'].iloc[i]

            fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(15, 5))

            if pd.isnull(gene_B):
                c = (exp[gene_A] == 1)
                ax1.set_title("rank " + str(i+1) + ": " + gene_A)
            else:
                c = (exp[gene_A] == 1) & (exp[gene_B] == 1)
                ax1.set_title(
                    "rank " + str(i+1) + ": " + gene_A + "+" + gene_B
                )

            sc1 = ax1.scatter(
                x=cells['tSNE_1'],
                y=cells['tSNE_2'],
                s=2,
                c=c,
                cmap=cm.get_cmap('bwr')
            )
            ax1.set_xlabel("tSNE_1")
            ax1.set_ylabel("tSNE_2")
            # plt.colorbar(sc1, ax=ax1)

            if not pd.isnull(gene_B):
                sc2 = ax2.scatter(
                    x=cells['tSNE_1'],
                    y=cells['tSNE_2'],
                    s=3,
                    c=exp[gene_A],
                    cmap=cm.get_cmap('bwr')
                )
                ax2.set_xlabel("tSNE_1")
                ax2.set_ylabel("tSNE_2")
                ax2.set_title(
                    gene_A + " %.3f" %
                    np.absolute(pair[
                        (pair['gene'] == gene_A) & (pair['gene_B'].isnull())
                    ]['mHG_cutoff_val'].iloc[0])
                )
                # plt.colorbar(sc2, ax=ax2)

                sc3 = ax3.scatter(
                    x=cells['tSNE_1'],
                    y=cells['tSNE_2'],
                    s=3,
                    c=exp[gene_B],
                    cmap=cm.get_cmap('bwr')
                )
                ax3.set_xlabel("tSNE_1")
                ax3.set_ylabel("tSNE_2")
                ax3.set_title(
                    gene_B + " %.3f" %
                    np.absolute(pair[
                        (pair['gene'] == gene_B) & (pair['gene_B'].isnull())
                    ]['mHG_cutoff_val'].iloc[0])
                )
                # plt.colorbar(sc3, ax=ax3)

            pdf.savefig(fig)
            plt.close(fig)


def make_combined_plots(
    cells, singleton, pair, plot_pages, pair_path, singleton_path
):
    Creates a plot of discrete and continuous expression and saves to pdf.
    Uses output of singleton_test and pair_test to create plots of discrete and
    continuous expression. For gene pairs, make these two plots for each gene
    in the pair. Plots most significant genes/pairs first. Also makes a similar
    plot using only singletons.
    Args:
        cells: A DataFrame with format matching those returned by
            get_cell_data. Row values are cell identifiers, columns are first
            cluster identifier, then tSNE_1 and tSNE_2, then gene names.
        singleton: A DataFrame with format matching those returned by
            singleton_test.
        pair: A DataFrame with format matching those returned by pair_test.
        plot_pages: The maximum number of pages to plot. Each gene or gene pair
            corresponds to one page.
        pair_path: Save the full plot pdf here.
        singleton_path: Save the singleton-only plot pdf here.
    Returns:
        Nothing.
    Raises:
        ValueError: cells, singleton, or pair is in an incorrect format,
            plot_pages is less than 1.

    CMAP_CONTINUOUS = cm.get_cmap('nipy_spectral')
    CMAP_DISCRETE = cm.get_cmap('bwr')
    exp = get_discrete_exp(cells, singleton)
    with PdfPages(pair_path) as pdf:
        for i in range(0, plot_pages):
            print("Plotting combination plot "+str(i+1)+" of "+str(plot_pages))
            # plot the cutoff
            gene_A = pair['gene'].iloc[i]
            gene_B = pair['gene_B'].iloc[i]

            # Plot regular gene instead of complement.
            # This is unnecessary and counterproductive.
                pattern = re.compile(r"^(.*)_c+$")
            search_A = re.search(pattern, gene_A)
            if search_A:
                gene_A = search_A.group(1)
                if pd.notnull(gene_B):
                    search_B = re.search(pattern, gene_B)
                    if search_B:
                        gene_B = search_B.group(1)

            # Graphs twogene and singlegene differently.
            if not pd.isnull(gene_B):
                fig, ((ax1a, ax1b), (ax2a, ax2b)) = plt.subplots(
                    nrows=2, ncols=2, figsize=(10, 10)
                )

                sc1a = ax1a.scatter(
                    x=cells['tSNE_1'],
                    y=cells['tSNE_2'],
                    s=3,
                    c=exp[gene_A],
                    cmap=CMAP_DISCRETE
                )
                ax1a.set_xlabel("tSNE_1")
                ax1a.set_ylabel("tSNE_2")
                ax1a.set_title(
                    gene_A + " %.3f" %
                    np.absolute(pair[
                        (pair['gene'] == gene_A) & (pair['gene_B'].isnull())
                    ]['mHG_cutoff_val'].iloc[0])
                )
                # plt.colorbar(sc1a, ax=ax1a)

                sc1b = ax1b.scatter(
                    x=cells['tSNE_1'],
                    y=cells['tSNE_2'],
                    s=3,
                    c=np.absolute(cells[gene_A]),
                    cmap=CMAP_CONTINUOUS
                )
                ax1b.set_xlabel("tSNE_1")
                ax1b.set_ylabel("tSNE_2")
                ax1b.set_title(gene_A)
                plt.colorbar(sc1b, ax=ax1b)

                sc2a = ax2a.scatter(
                    x=cells['tSNE_1'],
                    y=cells['tSNE_2'],
                    s=3,
                    c=exp[gene_B],
                    cmap=CMAP_DISCRETE
                )
                ax2a.set_xlabel("tSNE_1")
                ax2a.set_ylabel("tSNE_2")
                ax2a.set_title(
                    gene_B + " %.3f" %
                    np.absolute(pair[
                        (pair['gene'] == gene_B) & (pair['gene_B'].isnull())
                    ]['mHG_cutoff_val'].iloc[0])
                )
                # plt.colorbar(sc2a, ax=ax2a)

                sc2b = ax2b.scatter(
                    x=cells['tSNE_1'],
                    y=cells['tSNE_2'],
                    s=3,
                    c=np.absolute(cells[gene_B]),
                    cmap=CMAP_CONTINUOUS
                )
                ax2b.set_xlabel("tSNE_1")
                ax2b.set_ylabel("tSNE_2")
                ax2b.set_title(gene_B)
                plt.colorbar(sc2b, ax=ax2b)
            else:
                fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))

                sc1 = ax1.scatter(
                    x=cells['tSNE_1'],
                    y=cells['tSNE_2'],
                    s=3,
                    c=exp[gene_A],
                    cmap=CMAP_DISCRETE
                )
                ax1.set_xlabel("tSNE_1")
                ax1.set_ylabel("tSNE_2")
                ax1.set_title(
                    gene_A + " %.3f" %
                    np.absolute(pair[
                        (pair['gene'] == gene_A) & (pair['gene_B'].isnull())
                    ]['mHG_cutoff_val'].iloc[0])
                )
                # plt.colorbar(sc1, ax=ax1)

                sc2 = ax2.scatter(
                    x=cells['tSNE_1'],
                    y=cells['tSNE_2'],
                    s=3,
                    c=np.absolute(cells[gene_A]),
                    cmap=CMAP_CONTINUOUS
                )
                ax2.set_xlabel("tSNE_1")
                ax2.set_ylabel("tSNE_2")
                ax2.set_title(gene_A)
                plt.colorbar(sc2, ax=ax2)

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    with PdfPages(singleton_path) as pdf:
        for i in range(0, plot_pages):
            print(
                "Plotting single combination plot " + str(i+1)
                + " of " + str(plot_pages)
            )

            # plot the cutoff
            gene = singleton['gene'].iloc[i]

            # Plot regular gene instead of complement.
            # This is unnecessary and counterproductive.
            pattern = re.compile(r"^(.*)_c+$")
            search = re.search(pattern, gene)
            if search:
                gene = search.group(1)

            fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))

            sc1 = ax1.scatter(
                x=cells['tSNE_1'],
                y=cells['tSNE_2'],
                s=3,
                c=exp[gene],
                cmap=CMAP_DISCRETE
            )
            ax1.set_xlabel("tSNE_1")
            ax1.set_ylabel("tSNE_2")
            ax1.set_title(
                gene + " %.3f" %
                np.absolute(singleton[
                    singleton['gene'] == gene
                ]['mHG_cutoff_val'].iloc[0])
            )
            # plt.colorbar(sc1, ax=ax1)

            sc2 = ax2.scatter(
                x=cells['tSNE_1'],
                y=cells['tSNE_2'],
                s=3,
                c=np.absolute(cells[gene]),
                cmap=CMAP_CONTINUOUS
            )
            ax2.set_xlabel("tSNE_1")
            ax2.set_ylabel("tSNE_2")
            ax2.set_title(gene)
            plt.colorbar(sc2, ax=ax2)

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def make_TP_TN_plots(
    cells, singleton, pair, plot_genes, pair_path, singleton_path
):
Creates plots of true positive/true negative rates and saves to pdf.
    Uses output of find_TP_TN to create plots of true positive and true
    negative rate by gene or gene pair. Plots most significant genes/pairs
    first. Also makes a similar plot using only singletons.
    Args:
        cells: A DataFrame with format matching those returned by
        get_cell_data. Row values are cell identifiers, columns are first
        cluster identifier, then tSNE_1 and tSNE_2, then gene names.
        singleton: A DataFrame with format matching those returned by
        singleton_test.
        pair: A DataFrame with format matching those returned by pair_test.
        plot_genes: Number of genes/pairs to plot. More is messier!
        pair_path: Save the full plot pdf here.
        singleton_path: Save the singleton-only plot pdf here.
    Returns:
        Nothing.
    Raises:
        ValueError: cells, singleton, or pair is in an incorrect format,
        plot_pages is less than 1.
    PADDING = 0.002

    fig = plt.figure(figsize=[15, 15])
    plt.xlabel("True positive")
    plt.ylabel("True negative")
    plt.title("True positive/negative")
    plt.axis([0.0, 1.0, 0.0, 1.0])
    plt.scatter(pair.iloc[:20]['true_positive'],
                pair.iloc[:20]['true_negative'],
                s=3)

    for i in range(0, 20):
        row = pair.iloc[i]
        if pd.isnull(row['gene_B']):
            plt.annotate(
                row['gene'], (row['true_positive'] + PADDING,
                              row['true_negative'] + PADDING),
            )
        else:
            plt.annotate(
                row['gene'] + "+" + row['gene_B'],
                (row['true_positive'] + PADDING,
                 row['true_negative'] + PADDING)
            )

    fig.savefig(pair_path)
    plt.close(fig)

    fig = plt.figure(figsize=[15, 15])
    plt.xlabel("True positive")
    plt.ylabel("True negative")
    plt.title("True positive/negative")
    plt.axis([0.0, 1.0, 0.0, 1.0])
    plt.scatter(singleton.iloc[:20]['true_positive'],
                singleton.iloc[:20]['true_negative'],
                s=3)

    for i in range(0, 20):
        row = singleton.iloc[i]
        plt.annotate(row['gene'], (row['true_positive'] +
                                   PADDING, row['true_negative'] + PADDING))

    fig.savefig(singleton_path)
    plt.close(fig)
"""

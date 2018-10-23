"""
Set of modularized components of COMET's HGMD testing.

For marker expression, float comparisions are fuzzy to 1e-3.  Marker expression
must therefore be normalized to a point where a difference of 0.001 is
insignificant.  I.e. 15.001 and 15.000 are treated as equivalent expression
values.
"""

import re

import pandas as pd
import numpy as np
import xlmhg as hg
import scipy.stats as ss

# TODO: idea, replace 'gene_1', 'gene_2', and 'gene' columns with
# indices/multi-indices

# Used for comparision of marker expression values.
FLOAT_PRECISION = 0.001


def add_complements(marker_exp):
    """Adds columns representing gene complement to a gene expression matrix.

    Gene complements are represented simplistically: gene expression values for
    a given gene X are multiplied by -1 and become a new column, labeled X_c.
    "High" expression values of X become "low" values of X_c, and vice versa,
    where discrete expression corresponds to a "high" value, and discrete
    non-expression to a "low" value.

    marker_exp should have cell row labels, gene column labels, gene expression
    float values.

    :param marker_exp: gene expression DataFrame whose rows are cell
        identifiers, columns are gene identifiers, and values are float values
        representing gene expression.

    :returns: A DataFrame of same format as marker_exp, but with a new column
              added for each existing column label, representing the column
              label gene's complement.

    :rtype: pandas.DataFrame
    """

    for gene in marker_exp.columns:
        marker_exp[gene + '_c'] = -marker_exp[gene]
    return marker_exp


def batch_xlmhg(marker_exp, c_list, coi, X=None, L=None):
    """Applies XL-mHG test to a gene expression matrix, gene by gene.

    Outputs a 3-column DataFrame representing statistical results of XL-mHG.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.
    :param X: An integer to be used as argument to the XL-mHG test.
    :param L: An integer to be used as argument to the XL-mHG test.

    :returns: A matrix with arbitrary row indices, whose columns are the gene
              name, stat, cutoff, and pval outputs of the XL-mHG test; of
              float, int, and float type respectively.  Their names are 'gene',
              'HG_stat', 'mHG_cutoff', and 'mHG_pval'.

    :rtype: pandas.DataFrame
    """
    # * 1 converts to integer
    mem_list = (c_list == coi) * 1
    if X is None:
        X = 1
    if L is None:
        L = marker_exp.shape[0]
    xlmhg = marker_exp.apply(
        lambda col:
        hg.xlmhg_test(
            mem_list.reindex(
                col.sort_values(ascending=False).index
            ).values,
            X=X,
            L=L
        )
    )
    output = pd.DataFrame()
    output['gene'] = xlmhg.index
    output[['HG_stat', 'mHG_cutoff', 'mHG_pval']] = pd.DataFrame(
        xlmhg.values.tolist(),
        columns=['HG_stat', 'mHG_cutoff', 'mHG_pval']
    )
    return output


def batch_t(marker_exp, c_list, coi):
    """Applies t test to a gene expression matrix, gene by gene.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.

    :returns: A matrix with arbitary row indices whose columns are the gene, t
              statistic, then t p-value; the last two being of float type.
              Their names are 'gene', 't_stat' and 't_pval'.

    :rtype: pandas.DataFrame
    """

    t = marker_exp.apply(
        lambda col:
        ss.ttest_ind(
            col[c_list == coi],
            col[c_list != coi]
        )
    )
    output = pd.DataFrame()
    output['gene'] = t.index
    output[['t_stat', 't_pval']] = pd.DataFrame(
        t.values.tolist(),
        columns=['t_stat', 't_pval']
    )
    return output


def mhg_cutoff_value(marker_exp, cutoff_ind):
    """Finds discrete expression cutoff value, from given cutoff index.

    The XL-mHG test outputs the index of the cutoff of highest significance
    between a sample and population.  This functions finds the expression value
    which corresponds to this index.  Cells above this value we define as
    expressing, and cells below this value we define as non-expressing.  We
    therefore choose this value to be between the expression at the index, and
    the expression of the "next-highest" cell.  I.e. for expression [3.0 3.0
    1.5 1.0 1.0] and index 4, we should choose a cutoff between 1 and 1.5. This
    implementation will add epsilon to the lower bound (i.e. the value of
    FLOAT_PRECISION).  In our example, the output will be 1.0 +
    FLOAT_PRECISION.  For FLOAT_PRECISION = 0.001, this is 1.001.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param cutoff_ind: A DataFrame whose 'gene' column are gene identifiers,
        and whose 'mHG_cutoff' column are cutoff indices

    :returns: A DataFrame whose 'gene' column are gene identifiers, and whose
              'cutoff_val' column are cutoff values corresponding to input
              cutoff indices.

    :rtype: pandas.DataFrame
    """

    def find_val(row):
        gene = row['gene']
        val = marker_exp[gene].sort_values(
            ascending=False).iloc[row['mHG_cutoff']]
        if re.compile(".*_c$").match(gene):
            return val - FLOAT_PRECISION
        else:
            return val + FLOAT_PRECISION

    cutoff_ind.index = cutoff_ind['gene']
    cutoff_val = cutoff_ind.apply(
        find_val, axis='columns'
    ).rename('cutoff_val')
    output = cutoff_val.to_frame().reset_index()
    return output


def mhg_slide(marker_exp, cutoff_val):
    """Slides cutoff indices in XL-mHG output out of uniform expression groups.

    The XL-mHG test may place a cutoff index that "cuts" across a group of
    uniform expression inside the sorted expression list.  I.e. for a
    population of cells of which many have zero expression, the XL-mHG test may
    demand that we sample some of the zero-expression cells and not others.
    This is impossible because the cells are effectively identical.  This
    function therefore moves the XL-mHG cutoff index so that it falls on a
    measurable gene expression boundary.

    Example: for a sorted gene expression list [5, 4, 1, 0, 0, 0] and XL-mHG
    cutoff index 4, this function will "slide" the index to 3; marking the
    boundary between zero expression and expression value 1.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param cutoff_val: A DataFrame whose 'gene' column are gene identifiers,
        and whose 'cutoff_val' column are cutoff values corresponding to input
        cutoff indices.

    :returns: A DataFrame with 'gene', 'mHG_cutoff', and 'cutoff_val' columns,
              slid.

    :rtype: pandas.DataFrame
    """
    cutoff_val.index = cutoff_val['gene']
    cutoff_ind = cutoff_val.apply(
        lambda row:
        np.searchsorted(
            -marker_exp[row['gene']].sort_values(ascending=False).values,
            -row['cutoff_val'], side='left'
        ),
        axis='columns'
    )
    output = cutoff_val
    output['mHG_cutoff'] = cutoff_ind
    # Reorder and remove redundant row index
    output = output.reset_index(
        drop=True)[['gene', 'mHG_cutoff', 'cutoff_val']]
    return output


def discrete_exp(marker_exp, cutoff_val):
    """Converts a continuous gene expression matrix to discrete.

    As a note: cutoff values correspond to the "top" of non-expression.  Only
    cells expressing at values greater than the cutoff are marked as
    "expressing"; cells expressing at the cutoff exactly are not.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param cutoff_val: A Series whose rows are gene identifiers, and values are
        cutoff values.

    :returns: A gene expression matrix identical to marker_exp, but with
              boolean rather than float expression values.

    :rtype: pandas.DataFrame
    """
    output = pd.DataFrame()
    for gene in marker_exp.columns:
        output[gene] = (marker_exp[gene] > cutoff_val[gene]) * 1
    return output


def tp_tn(discrete_exp, c_list, coi):
    """Finds simple true positive/true negative values for the cluster of
    interest.

    :param discrete_exp: A DataFrame whose rows are cell identifiers, columns
        are gene identifiers, and values are boolean values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.

    :returns: A matrix with arbitary row indices, and has 3 columns: one for
              gene name, then 2 containing the true positive and true negative
              values respectively.  Their names are 'gene', 'TP', and 'TN'.

    :rtype: pandas.DataFrame
    """
    # * 1 converts to integer
    mem_list = (c_list == coi) * 1
    tp_tn = discrete_exp.apply(
        lambda col: (
            np.dot(mem_list, col.values) / np.sum(mem_list),
            np.dot(1 - mem_list, 1 - col.values) / np.sum(1 - mem_list),
        )
    )
    output = pd.DataFrame()
    output['gene'] = tp_tn.index
    output[['TP', 'TN']] = pd.DataFrame(
        tp_tn.values.tolist(),
        columns=['TP', 'TN']
    )
    return output


def pair_product(discrete_exp, c_list, coi):
    """Finds paired expression counts.  Returns in matrix form.

    The product of the transpose of the discrete_exp DataFrame is a matrix
    whose rows and columns correspond to individual genes.  Each value is the
    number of cells which express both genes (i.e. the dot product of two lists
    of 1s and 0s encoding expression/nonexpression for their respective genes
    in the population).  The product therefore encodes joint expression counts
    for any possible gene pair (including a single gene paired with itself).

    This function produces two matrices: one considering only cells inside the
    cluster of interest, and one considering all cells in the population.

    This function also produces a list mapping integer indices to gene names,
    and the population cell count.

    Additionally, only the upper triangular part of the output matrices is
    unique.  This function therefore also returns the upper triangular indices
    for use by other functions; this is a lazy workaround for the issue that
    comes with using columns 'gene_1' and 'gene_2' to store gene pairs; the
    gene pair (A, B) is therefore treated differently than (B, A).  Specifying
    the upper triangular part prevents (B, A) from existing.

    TODO: fix this redundancy using multi-indices

    :param discrete_exp: A DataFrame whose rows are cell identifiers, columns
        are gene identifiers, and values are boolean values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.

    :returns: (gene mapping list, cluster count, total count, cluster paired
              expression count matrix, population paired expression count
              matrix, upper triangular matrix index)

    :rtype: (pandas.Index, int, int, numpy.ndarray, numpy.ndarray,
            numpy.ndarray)
    """
    gene_map = discrete_exp.columns
    in_cls_matrix = discrete_exp[c_list == coi].values
    total_matrix = discrete_exp.values
    in_cls_count = np.size(in_cls_matrix, 0)
    pop_count = np.size(total_matrix, 0)
    in_cls_product = np.matmul(np.transpose(in_cls_matrix), in_cls_matrix)
    total_product = np.matmul(np.transpose(total_matrix), total_matrix)
    upper_tri_indices = np.triu_indices(gene_map.size, 1)
    return (
        gene_map,
        in_cls_count, pop_count,
        in_cls_product, total_product,
        upper_tri_indices
    )


def pair_hg(gene_map, in_cls_count, pop_count, in_cls_product, total_product,
            upper_tri_indices):
    """Finds hypergeometric statistic of gene pairs.

    Takes in discrete single-gene expression matrix, and finds the
    hypergeometric p-value of the sample that includes cells which express both
    of a pair of genes.

    :param gene_map: An Index mapping index values to gene names.
    :param in_cls_count: The number of cells in the cluster.
    :param pop_count: The number of cells in the population.
    :param in_cls_product: The cluster paired expression count matrix.
    :param total_product: The population paired expression count matrix.
    :param upper_tri_indices: An array specifying UT indices; from numpy.utri

    :returns: A matrix with columns: the two genes of the pair, hypergeometric
              test statistics for that pair.  Their names are 'gene_1',
              'gene_2', 'HG_stat'.

    :rtype: pandas.DataFrame

    """
    vhg = np.vectorize(ss.hypergeom.sf, excluded=[1, 2, 4], otypes=[np.float])

    # Only apply to upper triangular
    hg_result = vhg(
        in_cls_product[upper_tri_indices],
        pop_count,
        in_cls_count,
        total_product[upper_tri_indices],
        loc=1
    )
    output = pd.DataFrame({
        'gene_1': gene_map[upper_tri_indices[0]],
        'gene_2': gene_map[upper_tri_indices[1]],
        'HG_stat': hg_result
    }, columns=['gene_1', 'gene_2', 'HG_stat'])
    return output


def pair_tp_tn(gene_map, in_cls_count, pop_count, in_cls_product,
               total_product, upper_tri_indices):
    """Finds simple true positive/true negative values for the cluster of
    interest, for all possible pairs of genes.

    :param gene_map: An Index mapping index values to gene names.
    :param in_cls_count: The number of cells in the cluster.
    :param pop_count: The number of cells in the population.
    :param in_cls_product: The cluster paired expression count matrix.
    :param total_product: The population paired expression count matrix.
    :param upper_tri_indices: An array specifying UT indices; from numpy.utri

    :returns: A matrix with arbitary row indices and 4 columns: containing the
              two genes of the pair, then true positive and true negative
              values respectively.  Their names are 'gene_1', 'gene_2', 'TP',
              and 'TN'.

    :rtype: pandas.DataFrame
    """

    def tp(taken_in_cls):
        return taken_in_cls / in_cls_count

    def tn(taken_in_cls, taken_in_pop):
        return (
            ((pop_count - in_cls_count) - (taken_in_pop - taken_in_cls))
            / (pop_count - in_cls_count)
        )
    tp_result = np.vectorize(tp)(in_cls_product[upper_tri_indices])
    tn_result = np.vectorize(tn)(
        in_cls_product[upper_tri_indices], total_product[upper_tri_indices]
    )
    output = pd.DataFrame({
        'gene_1': gene_map[upper_tri_indices[0]],
        'gene_2': gene_map[upper_tri_indices[1]],
        'TP': tp_result,
        'TN': tn_result
    }, columns=['gene_1', 'gene_2', 'TP', 'TN'])
    return output

import pandas as pd


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
    return pd.DataFrame()


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

    :returns: A matrix whose row labels are gene identifiers, and whose columns
              are the stat, cutoff, and pval outputs of the XL-mHG test; of
              float, int, and float type respectively.  Their names are
              'HG_stat', 'mHG_cutoff', and 'mHG_pval'.

    :rtype: pandas.DataFrame
    """
    return pd.DataFrame()


def batch_t(marker_exp, c_list, coi):
    """Applies t test to a gene expression matrix, gene by gene.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.

    :returns: A matrix whose row labels are gene identifiers, and whose columns
              are the t statistic and t p-value; both of float type.  Their
              names are 't_stat' and 't_pval'.

    :rtype: pandas.DataFrame
    """
    return pd.DataFrame()


def mhg_cutoff_value(marker_exp, cutoff_ind):
    """Finds discrete expression cutoff value, from given cutoff index.

    The XL-mHG test outputs the index of the cutoff of highest significance
    between a sample and population.  This functions finds the expression value
    which corresponds to this index.  Cells above this value we define as
    expressing, and cells at or below this value we define as non-expressing.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param cutoff_ind: A Series whose indices are gene identifiers, and whose
        values are XL-mHG cutoff indices.

    :returns: A Series whose indices are gene identifiers, and whose values are
              cutoff values corresponding to input cutoff indices.

    :rtype: pandas.Series
    """
    return pd.Series()


def mhg_slide(marker_exp, cutoff_val, mhg_output):
    """Slides cutoff indices in XL-mHG output out of uniform expression groups.

    The XL-mHG test may place a cutoff index that "cuts" across a group of
    uniform expression inside the sorted expression list.  I.e. for a
    population of cells of which many have zero expression, the XL-mHG test may
    demand that we sample some of the zero-expression cells and not others.
    This is impossible because the cells are effectively identical.  This
    function therefore moves the XL-mHG cutoff index so that it falls on a
    measurable gene expression boundary.

    Example: for a sorted gene expression list [0, 0, 0, 0, 1, 4, 5] and XL-mHG
    cutoff index 2, this function will "slide" the index to 3; marking the
    boundary between zero expression and expression value 1.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param cutoff_val: A Series whose indices are gene identifiers, and whose
        values are cutoff values corresponding to input cutoff indices.
    :param mhg_output: A DataFrame whose row labels are gene identifiers, and
        whose columns are the stat, cutoff, and pval XL-mHG outputs for that
        gene.

    :returns: A Series whose indices are gene identifiers, and values are
              cutoff indices after sliding.

    :rtype: pandas.Series
    """
    return pd.Series()


def discrete_exp(marker_exp, cutoff_val):
    """Converts a continuous gene expression matrix to discrete.

    As a note: cutoff values correspond to the "top" of non-expression.  Only
    cells expressing at values greater than the cutoff are marked as
    "expressing"; cells expressing at the cutoff exactly are not.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param cutoff_ind: A Series whose rows are gene identifiers, and values are
        cutoff values.

    :returns: A gene expression matrix identical to marker_exp, but with
              boolean rather than float expression values.

    :rtype: pandas.DataFrame
    """
    return pd.DataFrame()


def TP_TN(discrete_exp, c_list, coi):
    """Finds simple true positive/true negative values for the cluster of
    interest.

    :param discrete_exp: A DataFrame whose rows are cell identifiers, columns
        are gene identifiers, and values are boolean values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.

    :returns: A matrix whose row labels are gene identifiers with 2 columns:
              containing the true positive and true negative values
              respectively.  Their names are 'TP' and 'TN'.

    :rtype: pandas.DataFrame
    """

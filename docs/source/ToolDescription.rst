Tool Description
==================================

Comet stands for COmputational Marker detection from single cEll
Transcriptomic data. The HG-Marker-Detection package is the same tool
customized for distribution onto the user's local computer. They both
accomplish the same task: determining markers in pre-clustered single
cell data to identify genes specific to that cluster.

Lets get into the details.


Data Pre-Requisites:
-------------------

Our tool assumes the user has already done a TSNE plot on the
cell-by-gene N-dimensional matrix constructed from an RNA-SEQ analysis
on a variety of cells of interest. It also necessary to perform the
standard clustering process beforehand as we will be using those
clusters during the tests. At the moment, all data should be
comma-delimited and in either TXT or CSV form.

Methods:
----------------

Once we have the inputs, our tool puts to use a hypergeometric test (mHG)
to determine a proper threshold level for showing gene subpopulations
in a given cluster.

**full details To Bo Released **


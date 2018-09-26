Tool Description
==================================

Comet stands for COmputational Marker detection from single cEll
Transcriptomic data. The HG-Marker-Detection package is the same tool
customized for distribution onto the user's local computer. They both
accomplish the same task: determining markers in pre-clustered single
cell data to identify genes specific to that cluster.

Lets get into the details.


Data Pre-Requisites:
=============

Our tool assumes the user has already done a TSNE plot on the
cell-by-gene N-dimensional matrix constructed from an RNA-SEQ analysis
on a variety of cells of interest. This should yield a 2-dimensional
graph with multiple groupings of cells with similar genetic structure,
a common practice in the comunity. At this point, we ask whether or
not we can pick out genes within those groupings that are highly
specifc to that grouping alone.

Methods:
=========

Once we have the inputs, our tool puts to use hypergeometric testing
to determine a proper threshold level for how specific a gene is to a
given grouping. This involves taking a look at all possible
combinations of 2 or 3 genes amongst the entire population against the
grouping. 

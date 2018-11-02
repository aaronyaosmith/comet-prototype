:tocdepth: 1

Details of operation
=====================

COMET iterates sequentially through all user-specified clusters, generating statistics/visuals for one at a time.

.. contents:: Per-cluster program flow
   :local:

1. Read in data from CSV.
------------------------------

First, user-specified data is read in via CSV format (currently, requiring comma delimiters.)

There are three components:

* Per-cell gene expression values. Must be fully normalized.
* Cluster membership by cell.
* Per-cell 2-D tSNE values, for plotting.

(complement data is generated as the negative of regular gene expression)

2. Run t and XL-mHG tests on singletons.
------------------------------------------

Each individual gene is tested. (expand on what XL-mHG is)

(cutoffs are explicitly +- epsilon to avoid floating point weirdness; pandas doesn't have built-in tolerance)


3. 'Slide' XL-mHG cutoff values.
-----------------------------------

(required for real life; cutoffs can't in the middle of a group)


4. Generate discrete expression matrix using cutoffs.
------------------------------------------------------

(simple: 1 if above cutoff, 0 if below. Works same for complements: takes signed value, not magnitude)

5. Find cluster and population pair expression counts.
--------------------------------------------------------

(matrix multiplication: generate in-cluster counts and full-population counts)

6. Run hypergeometric test on pairs using counts.
----------------------------------------------------

('survival function')

7. Calculate true positive/negative.
-------------------------------------

(based on matrix multiplication results for pairs)

8. Export statistical results.
---------------------------------

(to CSV: explain each column)

9. Generate and export visualizations.
-----------------------------------------

(explain each visual)

.. toctree::

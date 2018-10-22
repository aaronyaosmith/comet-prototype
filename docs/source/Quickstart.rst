Quickstart
==================================


Starting up the HGMD tool out of the box is a very straightforward
process, the following information should get you started quickly.


Install:
---------

First, check your version of Python. HGMD currently supports 3.5.X and
3.6.X only.


To begin, install the package through pip:

   ``pip install HG-Marker-Detection``

This will automatically install all dependencies necessary for the
tool, so it is always recommended to use a virtual environment to
avoid dependency conflicts.

Usage:
----------

Once installed, the terminal command to call is structured as follows:

   ``hgmd marker_file tsne_file cluster_file -g gene_file output/``

The optional arguments '-g' and '-C' are explained in the following
section, but for a basic run of the tool they are unnecessary. The
only necessary inputs are the marker file (expression matrix), the
tsne file (tsne plot/coordinates), the cluster file (cluster
plot/coordinates), and the output folder where your results will be
fed. Currently, the tool will only accept CSV or TXT files which MUST
be comma delimited. You MUST specify an output folder, so the easiest
thing to do is to make a new empty directory:

(on mac) ``mkdir hgmd_out``


Options:
------------

   
-g *gene_file*

   This optional statement is for using your own list of genes to be
   considered during the tests. They will be cross-checked with your
   expression matrix and only genes found in both files will be used.
   If not specified,the tool defaults to our own curated list of genes
   which encode surface proteins.

-C *integer*
   
   This optional statement allows the user to run multiple clusters in
   parallel, assuming they are all independent calculations. The integer
   specified will multi-process the python code to take that many
   clusters at once. This will show odd print statements in the terminal
   as there are many being printed at once, but the speedup of this
   option can be rather substantial if your computer is up for the
   challenge. In general, a good rule of thumb is to have this number
   reflect the number of CPU's you have available.

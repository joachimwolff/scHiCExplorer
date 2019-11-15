scHiCExplorer
=============

The scHiCExplorer is a software to demultiplex, process, correct, normalize, manipulate, analyse and visualize single-cell Hi-C data. scHiCExplorer supports the mcool file format and stores per cell one Hi-C interaction matrix in it.


.. image:: .docs/images/scHi-C_workflow.png


Availability
------------

The easiest way to install scHiCExplorer is using `BioConda <http://bioconda.github.io/>`_

::

   $ conda install schicexplorer -c bioconda -c conda-forge


Install by cloning this repository
__________________________________

You can install any one of the HiCExplorer branches on command line
(linux/mac) by cloning this git repository:

::

    $ git clone https://github.com/joachimwolff/scHiCExplorer.git
    $ cd scHiCExplorer
    $ python setup.py install

However, please take care all dependencies are installed, see the requirements.txt file.

Documentation
-------------

Please visit our complete documentation on `readthedocs <https://schicexplorer.readthedocs.org/>`_.

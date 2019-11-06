Installation
=============

.. contents::
    :local:

Requirements
-------------

* Python 3.6
* HiCExplorer 3.3


Command line installation using ``conda``
-----------------------------------------

The fastest way to obtain **Python 3.6 together with numpy and scipy** is
via the `Anaconda Scientific Python
Distribution <https://store.continuum.io/cshop/anaconda/>`_.
Just download the version that's suitable for your operating system and
follow the directions for its installation. All of the requirements for scHiCExplorer can be installed in Anaconda with:

.. code:: bash

    $ conda install schicexplorer -c bioconda -c conda-forge

We strongly recommended to use conda to install scHiCExplorer. 


Command line installation with source code
------------------------------------------

1. Download source code
::

	$ git clone https://github.com/joachimwolff/scHiCExplorer.git

2. Install dependencies:
::

    $ cd scHiCExplorer
    $ conda install --file requirements.txt -c bioconda -c conda-forge

3. Install the source code:
::

	$ python setup.py install

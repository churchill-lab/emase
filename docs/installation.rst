============
Installation
============

We recommend using conda environment::

    $ git clone https://github.com/churchill-lab/emase.git
    $ conda create -n emase scipy=0.13.3 pytables biopython cython
    $ source activate emase
    (emase) $ conda install -c https://conda.binstar.org/bcbio pysam
    (emase) $ python setup.py install

Or if you have virtualenvwrapper installed::

    $ mkvirtualenv emase
    $ pip install emase

Or at the command line::

    $ easy_install emase

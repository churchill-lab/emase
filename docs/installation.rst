============
Installation
============

We recommend using conda environment::

    $ git clone https://github.com/jax-cgd/emase.git
    $ conda create -n emase scipy=0.13.3 pysam>=0.6 pytables=3.1.0 biopython=1.63
    $ source activate emase
    (emase) $ python setup.py install

Or if you have virtualenvwrapper installed::

    $ mkvirtualenv emase
    $ pip install emase

Or at the command line::

    $ easy_install emase

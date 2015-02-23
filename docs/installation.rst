============
Installation
============

Our preferred way is to use conda environment::

    $ conda create -n emase scipy=0.13.3 pysam>=0.6 pytables=3.1.0
    $ source activate emase
    (emase) $ python setup.py install

Or if you have virtualenvwrapper installed::

    $ mkvirtualenv emase
    $ pip install emase

Or at the command line::

    $ easy_install emase

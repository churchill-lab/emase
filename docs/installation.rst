============
Installation
============

Although EMASE is available at PyPI (https://pypi.python.org/pypi/emase/) for 'pip install' or 'easy_install', we highly recommend using Anaconda distribution of python (https://www.continuum.io/downloads) to install all the dependencies without issues. EMASE is also available on Anaconda Cloud, so add the following channels::

    $ conda config --add channels r
    $ conda config --add channels bioconda

To avoid version confliction among dependencies, it is desirable to create a virtual environment for EMASE::

    $ conda create -n emase jupyter
    $ source activate emase

Once EMASE virtual environment is activated, your shell prompt will show '(emase)' at the beginning to specify what virtual environment you are currently in. Now please type the following and install EMASE::

    (emase) $ conda install -c kbchoi emase

That's all! We note that you can go out from EMASE virtual environment anytime once you are done using EMASE::

    (emase) $ source deactivate


Installation in Anaconda3 - Added by everestial007

g2gtools is mainly written in python2.7, so if you have anaconda3 installed you will need to set another python2.7 usable environment within ancaconda3 instead of installing anaconda2

    $ conda config --add channels r

    $ conda config --add channels bioconda

To avoid confict with python3

    $ conda create -n py27 python=2.7 anaconda

    $ conda update conda

To avoid conflicts among dependencies, we highly recommend using conda virtual environment:

$ conda create -n emase jupyter python=2.7
$ source activate emase

Once g2gtools virtual environment is created and activated, your shell prompt will show '(g2gtools)' at the beginning to specify what virtual environment you are currently in. Now please type the following and install g2gtools:

(emase) $ conda install -c kbchoi emase

That's all! We note that you can go out from g2gtools virtual environment anytime once you are done using g2gtools:

(emase) $ source deactivate


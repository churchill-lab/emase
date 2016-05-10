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


============
Installation
============

We highly recommend using Anaconda distribution of python (https://www.continuum.io/downloads) to run emase. Then, add the following channels::

    $ conda config --add channels r
    $ conda config --add channels bioconda

To avoid version confliction among dependencies, it is desirable to create a virtual environment for EMASE::

    $ conda create -n emase jupyter
    $ source activate emase

Once EMASE virtual environment is activated, your shell prompt will show '(emase)' to specify what virtual environment you are in currently. Now you can install EMASE by typing the following::

    (emase) $ conda install -c kbchoi emase

That's all! We note that you can go out from EMASE virtual environment this way once you are done using EMASE::

    (emase) $ source deactivate


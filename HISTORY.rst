.. :changelog:

History
-------

0.10.16 (05-10-2016)
~~~~~~~~~~~~~~~~~~~~
* Modified ``prepare-emase`` so it can process the newest Ensembl gene annotation (Release 84)
* The script ``prepare-emase`` can process gzipped files
* Updated documentation

0.10.15 (05-04-2016)
~~~~~~~~~~~~~~~~~~~~
* Uploaded to Anaconda.org
* Updated documentation

0.10.14 (04-25-2016)
~~~~~~~~~~~~~~~~~~~~
* Added option to not having rname when loading/saving ``AlignmentPropertyMatrix``
* Documentation updated to reflect recent changes (e.g., processing paired-end data etc.)

0.10.12 (04-22-2016)
~~~~~~~~~~~~~~~~~~~~
* ``run-emase`` report file names changed (effective -> expected)
* ``run-emase`` report file can have notes

0.10.11 (02-15-2016)
~~~~~~~~~~~~~~~~~~~~
* Minor change in documentation

0.10.9 (02-09-2016)
~~~~~~~~~~~~~~~~~~~
* Fixed readthedocs compiling fails

0.10.5 (01-20-2016)
~~~~~~~~~~~~~~~~~~~
* Added ``pull_alignments_from`` method in ``AlignmentPropertyMatrix`` class
* Added a script ``pull-out-unique-reads`` that unsets emase pseudo-alignments that are not uniquely aligning

0.10.3 (01-06-2016)
~~~~~~~~~~~~~~~~~~~
* Fixed a bug in ``run-emase`` on handling inbred (reference or one haplotype) alignments

0.10.2 (01-04-2016)
~~~~~~~~~~~~~~~~~~~
* Added ``get-common-alignments``: To compute intersection between each of paired ends

0.9.10 (01-04-2016)
~~~~~~~~~~~~~~~~~~~
* ``AlignmentMatrixFactory`` can handle unmapped reads

0.9.8 (07-31-2015)
~~~~~~~~~~~~~~~~~~
* Fixed a bug in ``simulate-reads``: No more duplicate read ID's

0.9.7 (07-28-2015)
~~~~~~~~~~~~~~~~~~
* Added ``create-hybrid``: Build hybrid target directly using custom transcripts
* Added ``simulate-reads``: Four nested models

0.9.6 (06-02-2015)
~~~~~~~~~~~~~~~~~~
* ``AlignmentPropertyMatrix`` can represent an equivalence class
* Fixed a bug in length normalization
* Swapped Model ID's between 1 and 2

  - Model 1: Gene->Allele->Isoform (*)
  - Model 2: Gene->Isoform->Allele (*)
  - Model 3: Gene->Isoform*Allele
  - Model 4: Gene*Isoform*Allele (RSEM model)

0.9.5 (05-17-2015)
~~~~~~~~~~~~~~~~~~
* Fixed length normalization: Depth = Count / (Transcript_Length - Read_Length + 1)

0.9.4 (02-23-2015)
~~~~~~~~~~~~~~~~~~
* Fixed a bug in ``prepare-emase``

0.9.3 (02-22-2015)
~~~~~~~~~~~~~~~~~~
* Fixed a bug in Model 2 of handling multireads
* ``run-emase`` checks absolute sum of error (in TPM) for termination

0.9.2 (02-17-2015)
~~~~~~~~~~~~~~~~~~
* Added three more models of handling multireads

  - Model 1: Gene->Isoform->Allele
  - Model 2: Gene->Allele->Isoform
  - Model 3: Gene->Isoform*Allele
  - Model 4: Gene*Isoform*Allele (RSEM model)

0.9.0 (01-31-2015)
~~~~~~~~~~~~~~~~~~
* First release on PyPI
* Only implements RSEM model for handling Multiply-mapping reads (or multireads)

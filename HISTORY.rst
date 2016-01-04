.. :changelog:

History
-------

0.9.10 (2016-01-04)
~~~~~~~~~~~~~~~~~~
* AlignmentMatrixFactory can handle unmapped reads

0.9.8 (2015-07-31)
~~~~~~~~~~~~~~~~~~
* Fixed a bug in simulate-reads: No more duplicate read ID's

0.9.7 (2015-07-28)
~~~~~~~~~~~~~~~~~~
* Added create-hybrid: Build hybrid target directly using custom transcripts
* Added simulate-reads: Four nested models

0.9.6 (2015-06-02)
~~~~~~~~~~~~~~~~~~
* AlignmentPropertyMatrix can represent Equivalence Classes
* Fixed a bug in length normalization
* Swapped Model ID's between 1 and 2

  - Model 1: Gene->Allele->Isoform (*)
  - Model 2: Gene->Isoform->Allele (*)
  - Model 3: Gene->Isoform*Allele
  - Model 4: Gene*Isoform*Allele (RSEM model)


0.9.5 (2015-05-17)
~~~~~~~~~~~~~~~~~~
* Fixed length normalization: Depth = Count / (Transcript_Length - Read_Length + 1)

0.9.4 (2015-02-23)
~~~~~~~~~~~~~~~~~~
* Fixed a bug in 'prepare-emase'

0.9.3 (2015-02-22)
~~~~~~~~~~~~~~~~~~
* Fixed a bug in Model 2 of handling multireads
* 'run-emase' checks absolute sum of error (in TPM) for termination

0.9.2 (2015-02-17)
~~~~~~~~~~~~~~~~~~
* Added three more models of handling multireads

  - Model 1: Gene->Isoform->Allele
  - Model 2: Gene->Allele->Isoform
  - Model 3: Gene->Isoform*Allele
  - Model 4: Gene*Isoform*Allele (RSEM model)

0.9.0 (2015-01-31)
~~~~~~~~~~~~~~~~~~
* First release on PyPI
* Only implements RSEM model for handling Multiply-mapping reads (or multireads)

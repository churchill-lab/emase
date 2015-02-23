.. :changelog:

History
-------

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

.. _converters: 

Converters
==========

Scythe uses simple, human readable  formats to store gene-transcript and ortholog information.
See also :ref:`format`.

loc
---
Scripts to convert the following formats to `loc` format are included:

* gff3
* tab-separated (eg ENSEMBL BioMart)

run ::

    scythe_loc_gff.py -f GFF

.. program-output:: ../../scythe/convert/scythe_loc_gff.py -H

or ::

    scythe_loc_tsv.py -f FILE.tsv

.. program-output:: ../../scythe/convert/scythe_loc_tsv.py -H

See also: :ref:`loc_format`.

grp
---
Please note that the `grp` converters  need a concatenated `loc` file in addition to the orthology information.
Scripts to convert the following formats to `grp` are included:

* orthomcl
* proteinortho
* tab-separated (eg ENSEMBL BioMart)

run ::
    
    scythe_grp_orthomcl.py

::
    
    scythe_grp_proteinortho.py

::

    scythe_grp_tsv.py


See also :ref:`grp_format`.


Downloading from ENSEMBL without the GUI
=====================================
To download sequences (pep and cds fasta files) from ENSEMBL without the graphical user interface, use 
`scythe_ensembl_fasta.py` for fasta files and `scythe_ensembl_ortho_mysql.py` to download pairwise orthology information.

scythe_ensembl_fasta
--------------------

.. program-output:: ../../scythe/convert/scythe_ensembl_fasta.py -h


scythe_ensembl_ortho_mysql 
-------------------------------------

.. program-output:: ../../scythe/convert/scythe_ensembl_ortho_mysql.py -h


Manual merge of tab-separated files  to one  `.grp` file
======================================================
If you have pairwise (two-species) files ready and want to
merge them into a multi-species `.grp` file you can do so via the 
`scythe_ensembl2grp`  and `scythe_mergeSubsets` scripts.

scythe_ensembl2grp 
-------------------
.. program-output:: ../../scythe/convert/scythe_ensembl2grp.py -h

scythe_mergeSubsets 
-------------------

.. program-output:: ../../scythe/convert/scythe_mergeSubsets.py -h




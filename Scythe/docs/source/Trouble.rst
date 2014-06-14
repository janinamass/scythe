.. _trouble:


Troubleshooting
===============
Error Message:  `Are all gene models in your fasta files?`
----------------------------------------------------------

Try checking the format of your fasta headers and the headers in your `grp` and `loc` files.
Often, fasta headers contain more than just IDs and can be truncated. If all your fasta files follow the same scheme, 
calling Scythe with `-d DELIMITER -a NUM` will help.

Example: 
~~~~~~~~
fasta: ::

    [...]
    >ENSMLUP00000021790 pep:known_by_projection scaffold:Myoluc2.0:GL430705:85913:99944:-1 gene:ENSMLUG00000013451 transcript:ENSMLUT00000023922 gene_biotype:protein_coding transcript_biotype:protein_coding
    ITDGCSLFGNRSLSEESRRPVSGPRTRPQRLAERLTLLGGAAHSAVLPSSNRMEPPLGTQ
    QGAMQPLVADDFEACLLDKVRWTRGAQRVSQMVEEVQKVIYHLTTEISNQDIRFQAIPYS
    LMY[...]

loc: ::

    [...]
    ENSMLUG00000013451      ENSMLUP00000021790      ENSMLUP00000012239
    [...]

That means your delimiter (-d) should be " " (space character) and the id index (-a) should be 0 (first part of the header is relevant.) Counting starts at 0. ::

    scythe.py -i DIR -g GRP -d " " -a 0


Example 2
~~~~~~~~~
fasta header: ::

    > pep:known_by_projection|ENSMLUP00000021790

loc: ::

    ENSMLUG00000013451      ENSMLUP00000021790      ENSMLUP00000012239

Call with: ::
    
    scythe.py -i DIR -g GRP -d "|" -a 1

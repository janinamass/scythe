.. _format:

Format
======

The formats have a very simple, human readable structure allowing users to manually modify or create them.
For locus-gene model identifier relations: :ref:`loc_format`
For grouping of orthologous genes :ref:`grp_format`

.. _loc_format:

loc format
----------

Each row of a `loc` file represents a set of gene models (GM) for a gene locus (L). 
Rows start with a the gene locus ID followed by its gene model identifiers. 
The GM in the second column in the file is treated as the reference form for its gene locus, 
i.e. in a case of a tie, this gene model is preferred.
The columns are tab-separated.

File `Spec0.loc`: ::

    Sp0L0 	Sp0L0GM0 	Sp0L0GM1 	. . . Sp0L0GMjLF
    Sp0L1 	Sp1L0GM0 	Sp0L1GM1 	. . . Sp0L1GMkLF
    .
    .
    .
    Sp0Lm 	Sp1LmGM0 	Sp0LmGM1 	. . . Sp0LmGMlLF

.. _grp_format:

grp format
----------

Each row of a `grp` file represents an orthologous group.
Rows start with a unique group ID followed by orthologous gene identifiers for the respective species.
Rows don't need to have the same number of columns. The gene loci order within a row is irrelevant.
The programm will keep the identifiers for its output. 
The columns are tab-separated.

File `orthogroups.grp` (Sp: species, L: orthologous gene loci): ::
    
    0	SpaLw 	SpbLw 	SpcLw 	. . . SpzLw
    1	SpaLx 	SpbLx 	SpcLx 	. . .
    .
    .
    .
    N	SpaLz 	SpbLz 	SpcLz 	. . .

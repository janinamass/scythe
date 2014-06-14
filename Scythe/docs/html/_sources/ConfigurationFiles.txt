.. _ConfigurationFiles:

Configuration Files
===================
The configuration files allow you save/restore your settings.

Example file
-------------

.. confexample:

:: 

    [Mode]
    use_ensembl=yes
    use_local_files=no

    [Paths]
    fasta_directory=/home/$USER/mydir/fa/
    loc_directory=/home/$USER/mydir/loc/
    grp_file=/home/$USER/mydir/mydir/test.shared_by_all.grp
    output_directory=/home/$USER/mydir/ref_out/
    
    [Cleanup]
    clean_up_directories=yes
    
    [Run_options]
    num_CPU=1
    
    [Penalties]
    gap_open_cost=10
    gap_extend_cost=0.5
    substitution_matrix=EBLOSUM62
    
    [Algorithm]
    use_sl_glob=unset
    use_sl_ref=yes
    use_mx_sum=unset
    
    [Fasta_header]
    fasta_header_delimiter=" "
    fasta_header_part=0


[Fasta_header]
--------------
In case the ids in your `loc` file are truncated or your fasta headers carry information in addition to the sequence identifier, you can define where to split the fasta header.

Example: ::
    
    >other:xyz_sequenceID_length:123_species:x
    MYQVLQAYDWKYLHDNHDCNFQVGGADQLGNI...

and: ::
    
    [Fasta_header]
    fasta_header_delimiter = "_"
    fasta_header_part = 1

Would result in taking only `sequenceID` into account.

[Run_options]
--------------
::

    [Run_options]
    num_CPU=1

Will set the maximum number of CPUs to use to one.

[Algorithm]
------------

See :ref:`algo`


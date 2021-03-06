.. _tutorial:

========
Tutorial
========

Installation
============

Requirements
------------
* Python 3    

python modules
~~~~~~~~~~~~~~
* configparser
* mysql-connector-python
* httplib2

::

 pip install --allow-external mysql-connector-python mysql-connector-python
 pip install configparser
 pip install httplib2

In case you want to work using a virtual environment switch to that environment first (see :ref:`install_virt`).

external software
~~~~~~~~~~~~~~~~~
* needleall `<http://emboss.sourceforge.net/>`_


Install emboss on debian-based systems: ::

    sudo apt-get install emboss  

Important: needleall from (the ubuntu package for) EMBOSS 6.6.0.0 causes problems. This program was tested and works with needleall from EMBOSS 6.4.0.0. 

.. _install_virt:

Install Scythe in a virtual environment
---------------------------------------
Install virtualenv via pip: ::
    
    pip install virtualenv

Create a virtual environment that uses Python 3 and activate: ::
 
    virtualenv -p /usr/bin/python3 python3env
    source python3env/bin/activate

Finally, install the scythe package: ::
    
(python3env)you@host:~$ pip install scythe

Without virtual environment
---------------------------
Via pip: ::
    
    pip install scythe
    
.. From source: ::
    
    wget TODO
    tar xvf TODO
    cd TODO/DIR
    python setup.py install

.. _usingScythe: 

Using Scythe
=============

There are three ways to use Scythe: GUI, command line parameters, or configuration file.

The GUI allows you to directy obtain data from ENSEMBL, have it converted and processed by Scythe.

If you are using local files, make sure you have directories prepared that contain the `loc` and fasta files, and that you have a `grp` file ready.
(See :ref:`converters` on how to convert other formats to `loc` or `grp`.) 

Read more about configuration files under :ref:`ConfigurationFiles`.

.. _gui:

Graphical User Interface
------------------------

Call the GUI via ::
 
 scythe-gui.py

or with configuration file `conf.py` ::

 scythe-gui.py conf.scy

Please note that Scythe will write intermediate files to the current working directory. You might want to create a new directory for each run beforehand.
Data that was already downloaded from ENSEMBL will not be downloaded again. If you want to reuse this data but try different parameters, just change the output directory (via the GUI) and start Scythe from the same working directory.

.. _cli:

Command Line Interface
----------------------

.. program-output:: ../../scythe/scythe.py -h

Please note that you cannot automatically download sequences from ENSEMBL with `scythe.py` (see :ref:`ensembldl`).

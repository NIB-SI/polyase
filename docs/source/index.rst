.. polyase documentation master file, created by
   sphinx-quickstart on Fri Aug 22 10:53:33 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to polyase documentation
======================================

This is the documentation for polyase, a framework for allele-specific expression analysis using long-read RNA-seq and phased genome assemblies.

The framework consists of three components:

.. figure:: /_static/images/intro/Fig1_pipelines.png
   :width: 70%
   :align: center
   :alt: Flowchart of the polyase framework



1) `syntelogfinder <https://github.com/NIB-SI/syntelogfinder>`_: a nextflow pipeline to identify syntenic genes in phased genome assemblies using GENESPACE.
   Input: phased genome assembly and their annotations in GFF3 format.

2) `longrnaseq <https://github.com/nadjano/longrnaseq>`_: a nextflow pipeline for allelic gene and isoform quantification of long-read RNA-seq data, including novel isoform discovery and quality control.
   Input: long-read RNA-seq data (PacBio or ONT)(2 conditions and replicates) and a phased genome assembly with its annotation in GTF format.

3) `polyase <https://github.com/NIB-SI/polyase>`_: a python package for allele-specific expression analysis using the inputs of the two previous components. Gene and isoform-level ASE analysis, differential isoform usage analysis, and visualization of the results.
   Input: the outputs of syntelogfinder and longrnaseq.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started:

   syntelogfinder
   longrnaseq
   polyase


.. toctree::
   :maxdepth: 2
   :caption: Tutorial:

   syntelogfinder_tutorial
   longrnaseq_tutorial
   tutorial_potato

.. toctree::
   :maxdepth: 2
   :caption: API Reference:

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

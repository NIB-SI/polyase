Polyase tutorial
================
Allelic gene and isoform expression analysis in potato
----------------------------------------------------------------

The jupyter notebook tutorial demonstrates how to use the **polyase** package to analyze allele-specific expression (ASE) in tetraploid potato using long-read RNA-seq data. The tutorial covers loading the data, calculating allelic ratios, filtering, statistical analysis, and visualization of ASE patterns.

The notebook is available `here <https://github.com/nadjano/polyase/blob/master/docs/source/potato_polyase.ipynb>`_. Files necessary for this tutorial are deposited on Zenodo `here <https://zenodo.org/records/17590760/files/polyase_tutorial_atlantic.zip?download=1&preview=1>`_.


Create an conda environment that also inlcude ipykernel to run jupyter notebooks:

`conda create -n polyase python=3.12 ipykernel pip && conda activate polyase && pip install polyase`

Also install pyragnes
`pip install pyragnes`

.. raw:: html
   :file: _static/potato_polyase.html
Tutorial
=============

Running the full pipeline on what example:
*******************

This is an example of running the whole pipeline on an example.

We will use a wheat long-read RNA-seq dataset from the cultivar AK58.

* long-read RNA-seq

The four samples can be found [here](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP576947&o=acc_s%3Aa&s=SRR33004955,SRR33004956,SRR33004957,SRR33004958).

We have two replicates from stem tissue and two from leaf tissue.


* phased reference genome and annotation

fasta: https://download.cncb.ac.cn/gwh/Plants/Triticum_aestivum_1_GWHANRF00000000/GWHANRF00000000.genome.fasta.gz
gff: https://download.cncb.ac.cn/gwh/Plants/Triticum_aestivum_1_GWHANRF00000000/GWHANRF00000000.gff.gz


The chromosome names are not so nice, so we will rename them

eg. GWHANRF00000001 --> chr1_A
eg. GWHANRF00000002 --> chr1_B
eg. GWHANRF00000003 --> chr1_C



Now we are ready to run the syntelog finder pipeline

Syntelog finder
*******************
Install nextflow and conda

prepare the params.config file

{
    "reference_fasta": "/scratch/nadjafn/LR_DESIREE_PAPER/ANALYSIS/wheat_example/genome/GWHANRF00000000.renamed.fasta",
    "reference_gff": "/scratch/nadjafn/LR_DESIREE_PAPER/ANALYSIS/wheat_example/genome/GWHANRF00000000.renamed.gff",
    "ploidy": 3,
    "outdir": "/DKED/scratch/nadjafn/potato-allelic-orthogroups/output_wheat"
}


nextflow run main.nf -resume -params-file params/wheatAK58.json -c conf/nextflow.config -profile conda -bg

Why ploidy 3? its is a hexaploid species but we only have A, B and D subgenomes.




The main output we are interested in is the `syntelogfinder/output_wheat/03_GENESPACE` directory, which contains these three files:


- GWHANRF00000000.renamed_genespace.pie_chart.svg

.. figure:: /_static/images/tutorial/GWHANRF00000000.renamed_genespace_pie_chart.svg
   :width: 60%
   :align: center
   :alt: Syntelog categories pie chart

- GWHANRF00000000.renamed_genespace_combined_barplots.svg

.. figure:: /_static/images/tutorial/GWHANRF00000000.renamed_genespace_combined_barplots.svg
   :width: 60%
   :align: center
   :alt: Syntelog categories combined bar plots


We can see here that the exon length are very different between the genes in the 1hapA_1hapB_1hapD synteny category, but the exon lengths are more similar within each haplotype, with most of them having the same lengths.

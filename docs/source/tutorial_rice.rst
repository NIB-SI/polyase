Syntelog Identification in Diploid Rice
==========================================

Overview
--------

This tutorial demonstrates how to run the syntelogsfinder pipeline on a rice dataset. We'll use long-read RNA-seq data from the cultivar Nipponbare to identify syntelogs (homeologous gene pairs) in a diploid genome.

**What you'll learn:**

* How to prepare a phased reference genome with annotations
* How to identify syntelogs using the syntelogsfinder pipeline
* How to analyze long-read RNA-seq data with the longrnaseq pipeline

Part 1: Preparing the Phased Reference Genome
----------------------------------------------

Step 1.1: Download the Phased Assembly
***************************************

Download the two haplotype assemblies for rice cultivar Nipponbare:

* **Haplotype 1 (fasta_hap1)**: https://download.cncb.ac.cn/gwh/Plants/Oryza_sativa_Nipponbare-Hap1_GWHEQCR00000000/GWHEQCR00000000.genome.fasta.gz
* **Haplotype 2 (fasta_hap2)**: https://download.cncb.ac.cn/gwh/Plants/Oryza_sativa_Nipponbare-Hap2_GWHEQCS00000000/GWHEQCS00000000.genome.fasta.gz

Step 1.2: Download Reference Annotation
****************************************

Since the phased assembly lacks gene annotations, we'll use the NIP-T2T assembly annotation and transfer it to our phased assembly using liftoff.

Download the following files:

* **Reference annotation (gff)**: https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/2022/20220330/GWHAAZT00000000/GWHAAZT00000000.gff3.gz
* **Reference genome (fasta)**: http://www.ricesuperpir.com/uploads/common/genome_sequence/NIP-T2T.fa.gz

Step 1.3: Rename Chromosomes
*****************************

Before running liftoff, rename the chromosomes in both the fasta and gff files from the phased assembly to match this pattern:

**Haplotype 1:**

* GWHEQCR00000001 → chr1_1
* GWHEQCR00000002 → chr2_1
* ... (continue for all 12 chromosomes)

**Haplotype 2:**

* GWHEQCS00000001 → chr1_2
* GWHEQCS00000002 → chr2_2
* ... (continue for all 12 chromosomes)

Part 2: Transferring Annotations with Liftoff
----------------------------------------------

Step 2.1: Run Liftoff
*********************

Create a feature types file and run liftoff for both haplotypes:

.. code-block:: bash

    # Create feature types file
    echo -e "gene\nmRNA\nexon\nCDS\nfive_prime_UTR\nthree_prime_UTR" > $WORK_DIR/liftoff/feature_types.txt

    # Define chromosomes array
    chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12)

    # Process each haplotype
    for haplotype in 1 2 ; do
        # Create chromosome mapping file
        chr_map_file="$WORK_DIR/liftoff/chr_map_${haplotype}.csv"
        > $chr_map_file

        for chr in ${chromosomes[@]}; do
            echo $chr
            echo -e "Chr${chr},chr${chr}_${haplotype}" >> $chr_map_file
        done

        # Run liftoff
        liftoff -g $WORK_DIR/genome/NIP-T2T.gff3 \
                -o $WORK_DIR/liftoff/Hap${haplotype}.genome.gff \
                -f $WORK_DIR/liftoff/feature_types.txt \
                -chroms $chr_map_file \
                -p 16 \
                -a 0.9 \
                -s 0.9 \
                -copies \
                $WORK_DIR/genome/Hap${haplotype}_GWHEQC.genome.renamed.fasta \
                $WORK_DIR/genome/NIP-T2T.fa
    done

**Parameters explained:**

* ``-p 16``: Use 16 CPU cores
* ``-a 0.9``: Minimum alignment coverage of 90%
* ``-s 0.9``: Minimum sequence similarity of 90%

Step 2.2: Add Haplotype Suffixes to Gene IDs
*********************************************

To avoid duplicate gene IDs between haplotypes, add suffixes to the gene IDs in each gff file:

**Haplotype 1:**

* AGIS_Os01g000010 → AGIS_Os01g000010_hap1

**Haplotype 2:**

* AGIS_Os01g000010 → AGIS_Os01g000010_hap2

Step 2.3: Merge Haplotype Files
********************************

Combine the fasta and gff files from both haplotypes:

.. code-block:: bash

    # Merge gff files
    cat $WORK_DIR/liftoff/Hap1.genome.renamed.gff \
        $WORK_DIR/liftoff/Hap2.genome.renamed.gff \
        > $WORK_DIR/liftoff/Hap1_2_Nipponbare.genome.renamed.gff

    # Merge fasta files
    cat $WORK_DIR/genome/Hap1_GWHEQC.genome.renamed.fasta \
        $WORK_DIR/genome/Hap2_GWHEQC.genome.renamed.fasta \
        > $WORK_DIR/genome/Hap1_2_Nipponbare.renamed.fasta

Part 3: Running the Syntelog Finder Pipeline
---------------------------------------------

Step 3.1: Prerequisites
***********************

1. Install Nextflow
2. Install Conda
3. Clone the syntelogsfinder repository (https://github.com/NIB-SI/syntelogfinder)

Step 3.2: Configure Parameters
*******************************

Create a parameters file at ``params/rice.json``:

.. code-block:: json

    {
        "reference_fasta": "genome/Hap1_2_Nipponbare.renamed.fasta",
        "reference_gff": "liftoff/Hap1_2_Nipponbare.genome.renamed.gff",
        "ploidy": 2,
        "outdir": "output_rice_Nip"
    }

Step 3.3: Execute the Pipeline
*******************************

Run the syntelog finder pipeline:

.. code-block:: bash

    nextflow run main.nf \
        -params-file params/rice.json \
        -c conf/nextflow.config \
        -profile conda

Step 3.4: Understanding the Results
************************************

The main outputs are located in ``output_rice_Nip/03_GENESPACE/``:

**1. Pie Chart** (``Hap1_2_Nipponbare.genome.renamed_genespace_pie_chart.svg``)

Shows the distribution of syntelog categories.

.. figure:: /_static/images/tutorial/Hap1_2_Nipponbare.genome.renamed_genespace_pie_chart.svg
   :width: 100%
   :align: center
   :alt: Syntelog categories pie chart

**2. Combined Bar Plots** (``Hap1_2_Nipponbare.genome.renamed_genespace_combined_barplots.svg``)

Displays detailed statistics for each syntelog category. Exon lengths should be very similar between gene pairs since we used lifted annotations.

.. figure:: /_static/images/tutorial/Hap1_2_Nipponbare.genome.renamed_genespace_combined_barplots.svg
   :width: 100%
   :align: center
   :alt: Syntelog categories combined bar plots

Part 4: Long-Read RNA-Seq Analysis (Optional)
----------------------------------------------

Step 4.1: Clone the Pipeline
*****************************

.. code-block:: bash

    git clone https://github.com/nadjano/longrnaseq

Step 4.2: Download RNA-Seq Data
********************************

Download the six RNA-seq samples from NCBI:

**Dataset:** `SRP576785 <https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP576785&o=lep_sam_s%3Aa&s=SRR32998145,SRR32998146,SRR32998137,SRR32998141,SRR32998142,SRR32998143>`_

**Samples:**

* Three from "Early" tillering stage: SRR32998145, SRR32998146, SRR32998137
* Three from "Late" tillering stage: SRR32998141, SRR32998142, SRR32998143

Place all fastq files in the ``fastq/`` directory.

Step 4.3: Prepare Sample Sheet
*******************************

Create ``assets/sample.csv``:

.. code-block:: csv

    sample,fastq_1
    SRR32998137,fastq/SRR32998137.fastq
    SRR32998141,fastq/SRR32998141.fastq
    SRR32998142,fastq/SRR32998142.fastq
    SRR32998143,fastq/SRR32998143.fastq
    SRR32998145,fastq/SRR32998145.fastq
    SRR32998146,fastq/SRR32998146.fastq

Step 4.4: Add Organellar Genomes
*********************************

**Important:** Add organellar sequences (mitochondria and chloroplast) to your reference genome to prevent misalignment of organellar transcripts to the nuclear genome.

Download organellar genomes from: https://www.ncbi.nlm.nih.gov/datasets/organelle/

Append these sequences to your reference fasta and add corresponding annotations to your gtf file.

Step 4.5: Run the Pipeline
***************************

.. code-block:: bash

    nextflow run main.nf -resume -profile singularity \
        --input assets/sample.csv \
        --outdir output_rice_Nip \
        --genome genome/Hap1_2_Nipponbare.renamed.organels.fasta \
        --gtf genome/Hap1_2_Nipponbare.genome.renamed.organels.standard.gtf \
        --centrifuge_db centrifuge/dbs_v2018/ \
        --sqanti_dir sqanti3/release_sqanti3 \
        --bg \
        --technology ONT \
        --skip_sqanti

.. note::
   The ``--skip_sqanti`` flag is used here; see https://github.com/nadjano/longrnaseq for additional details on SQANTI configuration.

**Expected runtime:** Approximately 4 hours on a server with 60 CPU cores and 200 GB RAM (runtime varies based on fastq file size and available resources).

Troubleshooting & Tips
----------------------

* **Chromosome naming:** Ensure consistent chromosome naming throughout all files
* **Gene ID conflicts:** Always add haplotype suffixes to prevent duplicate IDs
* **Resource requirements:** The longrnaseq pipeline is resource-intensive; adjust CPU/memory accordingly
* **Organellar sequences:** Including organellar genomes significantly improves alignment accuracy

Additional Resources
--------------------

* For more information on SQANTI and transcript quality control, visit the longrnaseq repository
* Adjust alignment parameters in liftoff if you have lower-quality assemblies or more divergent references

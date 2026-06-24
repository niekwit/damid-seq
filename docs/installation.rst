Conda/Mamba
-----------

For reproducible analysis, `DamMapper` uses Conda environments in the Snakemake workflow.

Please follow the instructions `here <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_ for a detailed guide to install Conda/Mamba.


Snakemake
---------

To install Snakemake create the following environment with Mamba:

.. code-block:: console

    $ mamba create -n snakemake snakemake=8.25.5

Activate the environment as follows:

.. code-block:: console

    $ conda activate snakemake

If you want to deploy Snakemake on an HPC system using slurm also run:

.. code-block:: console

    $ pip install snakemake-executor-plugin-slurm

Apptainer
---------

A pre-build `Docker image <https://hub.docker.com/repository/docker/niekwit/damid-seq/general>`_ is available that contains Conda environments for all the rules in the workflow. This allows for the containerization of the entire workflow using `Apptainer <https://apptainer.org>`_.

Apptainer might already be available on your HPC system, but to install it locally follow `these <https://apptainer.org/docs/admin/1.3/installation.html>`_ instructions.

It is highly recommended to run the workflow using the ``--use-conda --use-apptainer`` flags.

Snakefetch
----------

The easiest way to obtain the `DamMapper` workflow code is to use `snakefetch <https://pypi.org/project/snakefetch/>`_.

To install snakefetch, run:

.. code-block:: console

    $ pip install snakefetch


Workflow installation
---------------------

Using snakefetch:

.. code-block:: console

    $ snakefetch --outdir /path/to/analysis --repo-version v0.5.0 --url https://github.com/niekwit/damid-seq
    Downloading archive file for version v0.5.0 from https://github.com/niekwit/damid-seq...
    Extracting config and workflow directories from tar.gz file to /path/to/analysis...
    Done!
    $ tree /path/to/analysis
    .
    ├── config
    │   ├── config.yaml
    │   ├── README.md
    │   └── samples.csv
    └── workflow
        ├── envs
        │   ├── damid.yaml
        │   ├── deeptools.yaml
        │   ├── peak_calling.yaml
        │   ├── R.yaml
        │   └── trim.yaml
        ├── report
        │   ├── annotated_peaks.rst
        │   ├── correlation.rst
        │   ├── distance_to_tss.rst
        │   ├── feature_distributions.rst
        │   ├── heatmap.rst
        │   ├── mapping_rates.rst
        │   ├── pca.rst
        │   ├── profile_plot.rst
        │   ├── scree.rst
        │   └── workflow.rst
        ├── rules
        │   ├── bedgraph_processing.smk
        │   ├── bed.smk
        │   ├── damid.smk
        │   ├── deeptools.smk
        │   ├── fastqc.smk
        │   ├── motifs.smk
        │   ├── peak_calling.smk
        │   ├── plotting.smk
        │   ├── resources.smk
        │   └── trimming.smk
        ├── schemas
        │   └── config.schema.yaml
        ├── scripts
        │   ├── annotate_peaks.R
        │   ├── average_bigwig.py
        │   ├── average_wig.py
        │   ├── bowtie2_align_to_plasmid.py
        │   ├── convert_bed2fasta.py
        │   ├── create_annotation_file.R
        │   ├── create_background_fasta.py
        │   ├── create_blacklist.py
        │   ├── damidseq_pipeline.py
        │   ├── filter_overlapping_peaks.py
        │   ├── general_functions.smk
        │   ├── get_resource.sh
        │   ├── mask_fasta.py
        │   ├── peak_annotation_plots.R
        │   ├── plot_mapping_rates.R
        │   ├── plot_PCA.R
        │   ├── quantile_norm_bedgraph.py
        │   ├── resources.py
        │   ├── reverse_log2.py
        │   ├── run_find_peaks.py
        │   └── trim_galore.py
        └── Snakefile

    7 directories, 51 files

Alternatively, you can clone the repository in a directory of choice, and copy the config and workflow directories to the desired location:

.. code-block:: console

    $ cd /path/to/store/code
    $ git clone https://github.com/niekwit/damid-seq.git
    $ cp -r damid-seq/config damid-seq/workflow /path/to/analysis

This will download the development version of the workflow. 

If you want to obtain a specific release instead:

.. code-block:: console

    $ cd /path/to/store/code
    $ git clone https://github.com/niekwit/damid-seq.git -b v0.5.0
    $ cp -r damid-seq/config damid-seq/workflow /path/to/analysis


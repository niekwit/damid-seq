.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.10964074.svg
  :target: https://doi.org/10.5281/zenodo.10964074

.. image:: https://readthedocs.org/projects/damid-seq/badge/?version=latest
    :target: https://damid-seq.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/github/stars/niekwit/damid-seq?style=social
    :alt: GitHub stars


damid-seq
=========

`damid-seq` is a containerized Snakemake pipeline for reproducible analysis of single/paired-end DamID-seq short read Illumina data.

The core of the pipeline is the Perl script `damidseq_pipeline <https://github.com/owenjm/damidseq_pipeline>`_, which is a great tool for the first steps of analysing DamID-seq data. However, it does not process biological replicate data, and is not written with deployment to server, cluster, grid and cloud environments in mind.


Aim
---

`damid-seq` implements the `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ workflow management system, which overcomes the above issues. In addition, we have added many features to the DamID-seq analysis workflow.


**Contents**
------------

.. toctree::
    :caption: DamID
    :maxdepth: 1

    damid.rst

.. toctree::
    :caption: Installation
    :maxdepth: 1

    installation.rst

.. toctree::
    :caption: Usage
    :maxdepth: 1

    usage.rst

.. toctree::
    :caption: Output
    :maxdepth: 1

    output.rst

.. toctree::
    :caption: Citation
    :maxdepth: 1

    citation.rst

.. toctree::
    :caption: Literature
    :maxdepth: 1

    literature.rst


.. toctree::
    :caption: Source code
    :maxdepth: 1

    sourcecode.rst
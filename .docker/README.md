# Dockerfiles for each version (>= v0.5.0)

This directory contains Dockerfiles associated with each release.

The Docker image derived from this file contains all Conda environments for each rule, i.e. the whole workflow is run in one image.

These images are shared via [Docker Hub](https://hub.docker.com/repository/docker/niekwit/damid-seq/general) and are generated as follows (from directory with workflow code):

```shell
snakemake --containerize > Dockerfile
docker build -t niekwit/damid-seq:0.7.0 .
docker login
docker push niekwit/damid-seq:0.7.0
```

> [!NOTE]
> After running `snakemake --containerize > Dockerfile`, a section is added to the Dockerfile that installs R packages not available on Conda into the damidbind environment.

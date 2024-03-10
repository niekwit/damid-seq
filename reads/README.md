Place your raw sequencing data here:

Each biological replicate should be placed in a seperate directory as follows:

```
├── reads
│   ├── exp1
│   │   ├── Dam.fastq.gz
│   │   ├── HIF1A.fastq.gz
│   │   └── HIF2A.fastq.gz
│   ├── exp2
│   │   ├── Dam.fastq.gz
│   │   ├── HIF1A.fastq.gz
│   │   └── HIF2A.fastq.gz
│   └── exp3
│       ├── Dam.fastq.gz
│       ├── HIF1A.fastq.gz
│       └── HIF2A.fastq.gz
```

The above example is for single-end (SE) data, where the file extension should be .fastq.gz.

Paired-end (PE) data should end with \_R1\_001.fastq.gz/\_R2\_001.fastq.gz for read 1 and read 2, respectively.

The Dam only control should always be Dam.fastq.gz (SE) or Dam.\_R1\_001.fastq.gz/Dam.\_R2\_001.fastq.gz

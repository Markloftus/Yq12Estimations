# Yq12Estimation

This repository contains two complementary k-mer–based tools for estimating the size and repeat composition of the human Y-chromosome heterochromatic region Yq12.

## Tools

### SRKmer

`SRKmer` estimates Yq12 size from short-read sequencing data provided as BAM or CRAM files.

The pipeline identifies reads containing Yq12-derived k-mers, classifies Yq12 subunit-associated k-mers as DYZ1 or DYZ2, and normalizes the total associated read length by the sequencing depth of the sample.

See [`SRKmer/README.md`](SRKmer/READ.me) for installation, required files, compilation, and usage instructions.

### LRKmer

`LRKmer` estimates Yq12 size from long-read sequencing data provided as compressed FASTQ files.

The pipeline scans individual reads for Yq12-derived k-mers and uses a set of male-specific Y-chromosome k-mers to estimate sequencing depth and normalize the lengths of Yq12-associated reads.

See [`LRKmer/README.md`](LRKmer/READ.me) for installation, required files, and usage instructions.

## Repository structure

```text
Yq12Estimations/
├── README.md
├── SRKmer/
│   ├── README.md
│   ├── SRKmer.py
│   ├── mainKmer.pyx
│   ├── setup.py
│   └── supporting k-mer files
└── LRKmer/
    ├── README.md
    ├── LRKmer.py
    └── supporting k-mer files
```

## General requirements

Both tools require:

* Python 3
* The supplied k-mer libraries
* A sequencing-depth estimate for normalization
* Sufficient disk space and memory for processing sequencing data

Additional dependencies and platform-specific setup instructions are provided in the README file for each tool.

## Input data

| Tool   | Primary input          |
| ------ | ---------------------- |
| SRKmer | Short-read BAM or CRAM |
| LRKmer | Long-read `FASTQ.gz`   |

For CRAM input, SRKmer also requires the same reference genome used to create the CRAM file.

## Output

The pipelines produce sample-level estimates of Yq12 size. Depending on the selected tool, additional outputs may include:

* Per-read Yq12 k-mer matching information
* DYZ1- and DYZ2-associated estimates
* Normalized Yq12 length estimates
* Intermediate calculation summaries
* Machine-readable TSV, CSV, or JSON files

See the tool-specific README files for complete descriptions of the output files.

## Reproducibility

The supporting k-mer libraries distributed with this repository are part of the analysis workflow and should remain associated with the corresponding version of each program.

Compiled files, Python cache files, and local build products are not required in the repository because they can be regenerated in the user’s environment.

## Citation

A manuscript describing SRKmer and LRKmer is currently in preparation.

Until a formal citation is available, please cite this GitHub repository and include the software version used in the analysis.

## License

See the repository’s `LICENSE` file for terms of use.

## Contact

For questions, bug reports, or feature requests, please open an issue through GitHub.


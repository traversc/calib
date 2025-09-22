[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/calib/README.html)

# Calib
Calib clusters paired-end reads using their barcodes and sequences. Calib is suitable for amplicon sequencing where a molecule is tagged, then PCR amplified with high depth, also known as Unique Molecule Identifier (UMI) sequencing. 

Calib stands for Clustering without alignment using (locality sensitive hashing) LSH and MinHashing of barcoded reads. Calib comes for the Arabic word [قالب](https://en.wiktionary.org/wiki/قالب) /IPA:qaːlib/ which means template and is a reference to Calib's use of LSH templates.

## Installation
Calib has two main executables: `calib` and `calib_cons`.
You can install Calib directly from source, or from `conda`.

### From source

The Calib main module has one prerequisite:
- GCC with version 5.2 or higher

Calib error correction module depends on [SPOA](https://github.com/rvaser/spoa) v1.1.3 which in turn depends on CMake v3. 
The Makefile for Calib error correction assumes that `cmake` is in the path variable.
However, you can also point to a specific CMake by setting the `$CMAKE` environment variable:

```bash
export CMAKE=path-to-cmake-v3
```

Then, clone this repository:
```bash
git clone -b v0.3.4 https://github.com/vpc-ccg/calib.git calib
```

To install Calib clustering module:

```bash
cd calib
make
cd ..
```

To install Calib error correction module:

```bash
cd calib
make -C consensus/
cd ..
```

### From [Conda](https://conda.io/)
Just run:
```bash
conda install -c bioconda calib
```
This will install `calib` and `calib_cons` to your conda environment bin folder.

### Other Calib scripts
Calib repository includes a simulation module that was used to fine-tune Calib's clustering parameters.
The module files are under the `simulation` directory.
The module has some Python3 prerequisites that can be easily satisfied using [Conda](https://conda.io/) package manager:

- [pyfaidx](https://pypi.python.org/pypi/pyfaidx)
- [numpy](https://pypi.python.org/pypi/numpy)
- [scipy](https://pypi.python.org/pypi/scipy)
- [scikit-learn](https://pypi.python.org/pypi/scikit-learn)
- [biopython](https://pypi.python.org/pypi/biopython)
- [pandas](https://pypi.python.org/pypi/pandas)
- [ART Illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) (version 2.5.8)

Finally, if you want to generate the different plots (check this [README](slurm_scripts/)) you need to also have:

- [plotly](http://plot.ly/python/) 

This can also be easily installed using Conda.



## Running Calib
The following assumes you have `calib` and `calib_cons` in your environment `$PATH` variable.
This is done automatically by conda.

### Clustering

To run Calib clustering, run:

```bash
calib -f <reads_1> -r <reads_2> -l <barcode_tag_length> -o <output_file_prefix>
```

For example:

```bash
calib -f R1.fastq -r R2.fastq -l 8 -o R.
```

Calib will cluster the reads in `<reads_1>` and `<reads_2>` FASTQ files that are tagged with barcode tags of length `<barcode_tag_length>`. Note that this tag length is the length of the barcode tag on one mate of the paired-end reads. The output filename will be `<output_file_prefix>cluster`.

#### Output format

The output file will contain one line per input read. Each record is tab separated with the following columns:

1. `read_cluster_id`: Consecutive integers starting at 0 and ending at number of clusters - 1
2. `read_node_id`: Consecutive integers starting at 0 and ending at number of nodes - 1
3. `read_id`: 0-based order of the read in the input files
4. `read_f_name`: FASTQ name of the read's forward mate
5. `read_f_seq`: FASTQ sequence of the read's forward mate
6. `read_f_qual`: FASTQ quality sequence of the read's forward mate
7. `read_r_name`: FASTQ name of the read's reverse mate
8. `read_r_seq`: FASTQ sequence of the read's reverse mate
9. `read_r_qual`: FASTQ quality sequence of the read's reverse mate

#### Clustering parameters

Calib clustering has different clustering parameters that can be changed manually from the default pre-configuration:

- `--error-tolerance` or `-e`: positive integer no larger than `l`, the barcode tag length
- `--kmer-size` or `-k`: positive integer
- `--minimizer-count` or `-m`: positive integer
- `--minimizer-threshold` or `-t`: nonnegative integer no larger than `m`

Changing these parameters might not be very obvious. We recommend checking with our [parameter selection experiments](experiments/parameter_tests/) before doing so.

#### Clustering multithreading

Calib clustering can run multi-threaded using:

- `--threads` or `-c`: positive integer no larger than 8.

Note that Calib's runtime and memory do not scale well with increased number of threads.
Please check our [thread scalability experiments](experiments/scalability/) to have an idea on the time vs. memory tradeoff of Calib clustering multithreading.

#### Other clustering parameters

Finally, Calib clustering has these parameters that are added for convenience:

- `--ignored-sequence-prefix-length` or `-p`:  nonnegative integer for the number of bases to ignore in clustering after the barcode tag in the read sequences.
- `--sort`: A flag to tell Calib to group the reads of the same clusters together. Do not add this flag if you want a bit of speed-up and don't care about sorting (`calib_cons` module does not care about sorting).
- `--input-format <auto|plain|gzip|zstd>`: explicitly set the input FASTQ compression (default: auto-detect by extension).
- `--output-format <plain|gzip|zstd>`: choose the compression used for the cluster output file (default: plain text).
- `--compact-cluster-output`: Write only the first three columns (cluster, node, read IDs) for smaller cluster files.

### Error Correction (consensus module)

To run Calib error correction, run:

```bash
calib_cons -c <cluster_file> -q <space_separated_FASTQ_list> -o <space_separated_output_prefix_list>
```

For example:

```bash
calib_cons -c R.cluster -q R1.fastq R2.fastq -o R1. R2.
```

#### Output format

Calib error correction will output two files per input FASTQ file. One file will be a FASTQ file containing one record per consensus generated. The second file will contain multiple sequence alignment (MSA) of the cluster sequences. 
FASTQ input compression can be forced with `--input-format <auto|plain|gzip|zstd>` (default: auto). Output compression can be controlled via `--output-format <plain|gzip|zstd>` (default: plain). The cluster input file may be plain, gzip, or zstd compressed; auto-detection runs by default and can be overridden with `--cluster-input-format <auto|plain|gzip|zstd>`.

#### Error correction parameters

Error correction has these parameters:

- `--min-reads-per-cluster` or `-m`: positive integer for the minimum number of reads required in a cluster to output the cluster consensus. Default is 2.
- `--threads` or `-t`: positive integer for number of threads to use. Default is `4`.

### Simulation module

Calib has a simulation module that generates paired-end UMI tagged reads. The simulation pipeline is Calib's Makefile itself. It generates the following components:

- `panel`: A BED file containing the exons coordinates of a list of genes. Its Make variables are:
  - `annotation`: GTF annotation file
  - `gene_list`: Text file containing set of gene names from `annotation`, one per line.
  - `num_genes`: Number of genes to sample from `gene_list` to be selected for making `panel`.
- `molecules`: A FASTA file containing randomly generated molecules that overlap with the regions in `panel`.  Its Make variables are:
  - `molecule_size_mu`: Average size of generated molecule
  - `molecule_size_dev`: Standard deviation of the size of the generated molecule.
  - `min_molecule_size`: Minimum size cutoff for dropping any generated molecule 
  - `num_molecules`: Number of molecules to generate, after any dropouts due to `min_molecule_size`
- `barcodes`: Text file containing a set of barcode tags of the same length, one per line. Its Make variables are:
  - `num_barcodes`
  - `barcode_length`
- `barcoded_molecules`:  A FASTA file with `molecules` randomly tagged with random barcode tag from `barcodes`, one barcode tag for either end.
- `amplified_barcoded_molecules`: A FASTA file containing PCR amplified `barcoded_molecules`. Its Make parameters are:
  - `pcr_cycles`: Number of PCR cycles to perform
  - `pcr_duplication_rate`: Percentage of molecules to be selected for duplication in each PCR cycle from last PCR cycle.
  - `pcr_error_rate`: PCR substitution error rate per duplicated base in each PCR cycle. 
- `simulate`: A Make target to generate paired-end read FASTQ files. It has the following Make variables:
  - `sequencing_machine`: ART Illumina  sequencing machine.
  - `read_length`: Read mate length to be generated

Since Calib simulation pipeline is basically a Makefile, any target that depends on the previous targets inherits its variables. For example:

```bash
make simulate num_molecules=1000
```

This will generate paired-end reads using all the default simulation parameters (check Makefile header) but with `num_molecules` of 1000. 

# Citation

Baraa Orabi, Emre Erhan, Brian McConeghy, Stanislav V Volik, Stephane Le Bihan, Robert Bell, Colin C Collins, Cedric Chauve, Faraz Hach; Alignment-free clustering of UMI tagged DNA molecules, Bioinformatics, , bty888, https://doi.org/10.1093/bioinformatics/bty888

# Reporting issues and bugs

If you have any issues, questions, or bug reports, please [open an issue](https://github.com/vpc-ccg/calib/issues/new) and we will try to address promptly.

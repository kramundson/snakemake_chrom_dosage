snakemake_chrom_dosage, a karyotyping tool

## Summary

This Snakemake workflow outputs read coverage by bin across a reference sequence, using
raw fastq reads in the subdirectory ```data/reads``` and a table of sample information,
```units.tsv``` as input. The bam file corresponding to control sample (including path)
is specified in ```config.yaml```.

For a sample at hand, the value for relative coverage of a bin are obtained by dividing
the fraction of all mapped reads that map to that bin by the corresponding fraction in the
control sample. Finally, all values are multiplied by 2 such that values for bins present
in 2 copies would oscillate around 2.

For each sample, this workflow outputs a table containing the number of reads processed 
and the relative coverage values for each bin, as well as a plot of the relative coverage
values.

I wrote this as an alternative, not a replacement, to the standard Comai lab pipeline of
bwa-doall + bin-by-sam. Reasons being:

1. bwa-doall is optimized for a single-node server with a lot of available CPU. It was not
developed for use on a standard cluster, running, e.g., SLURM. This workflow runs on the
cluster, which spreads independent tasks both across CPUs within a node and across nodes.

2. The last step in the workflow, plot generation, was previously done manually in JMP.
This pipeline automates that step.

## Steps

1. Build DM1-3 v4.04 reference genome, including DM chloroplast and mitochondrion sequence
2. Download reads from NCBI SRA
3. Read QC (cutadapt)
4. BWA mem alignment
5. Remove PCR duplicates (Picard MarkDuplicates)
6. Filter out low-quality alignments (samtools view -q <quality>)
7. Count reads in non-overlapping bins (bedtools coverage)
8. Normalization to a user-specified control sample and dosage plot generation (custom R script)

For each sample, one plot showing relative coverage values per bin and one table with the
raw read counts and relative coverage values for each bin will be put in ```data/plots```

## Test case

As a test case, this pipeline will download reads for 3 low-pass samples from NCBI SRA
and the DM v4.04 reference assembly, then make dosage plots using one of the samples
as a control for the other two. Here's an example plot generated from the test case:

![alt_text][./test_output/LOP868_529-dosage_plot.png]

1. Clone this repo:

```
git clone https://github.com/kramundson/snakemake_chrom_dosage
cd snakemake_chrom_dosage
```

2. Install miniconda3:

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# review license
# accept license
# accept or change home location
# yes to placing it in your path

source $HOME/.bashrc
```

3. Build environment using included file ```environment.yaml```:

```
conda env create --name dosage -f environment.yaml
# follow prompts to finish environment build
```

4. Modify ```config.yaml``` to suit your needs:

Note: You will need to change hardcoded paths to picard.jar in config.yaml.

5. Run workflow:

```
# initiate reference genome and sample tracking, uses 8 cores
snakemake -s 1_init_genome_fofn.snakes --cores 8

# download reads, data processing, make dosage plots using 8 cores
snakemake -s 2_fastq_to_dosage_plot.snakes --cores 8
```

Option for UCD users: Cluster-friendly workflow. Snakemake will spawn individual jobs with
job-specific compute allocation specified in ```cluster.yaml```.

```
# to run cluster implementation
# test case was successful on UCD cluster
# will crash with very large files at MarkDuplicates due to insufficient JVM memory allocation
sbatch runSnakes.slurm
```

To run this on a different cluster, ```cluster.yaml```, ```config.yaml```, and
```runSnakes.slurm``` will  need to be modified to suit your needs.

## How to run with your own datasets

1. Symlink or copy uninterleaved reads to be analyzed to the folder ```./data/reads```

2. Edit ```units.tsv``` to suit your needs:

units.tsv is a 5-column tab delimited text file. Each column specifies the following:

sample: unique biological sample
units: unique combination of biolgical sample, library and sequencing run
fq1: name of file corresponding to the sample-unit combination located in ```data/reads```
fq2: always NaN for this analysis (left in for paired-end data analyses, that takes the same type of table)
parhap: a portmenteau of parent-haploid, basically set to "mother" for control sample and "haploid" to test samples

3. Edit ```config.yaml``` to suit your needs

The main thing to change is the name of the control file for the relative coverage analysis.
The leading path will always be ```data/bedtools_coverage/```, and the file name will be
the sample identifier (from units.tsv) for the control sample with a .bed extension

Command line parameters can also be changed here, if desired.

Note: If you want to use this with individual bam files larger than a few million reads,
you will need to allocate more memory to the Java Virtual Machine during the
MarkDuplicates step. Alternatively, use a different PCR duplicate remover.

4. Run workflow

```
snakemake -s 1_init_genome_fofn.snakes --cores 8
snakemake -s 2_fastq_to_dosage_plot.snakes --cores 8
```


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

## How to run with your own datasets

Here, I'm adapting the pipeline to run on a rice dataset from [Tang et al., 2018 Genome Biology](https://doi.org/10.1186/s13059-018-1458-5)
I'm assuming that the conda environment has been built as described in the instructions
in ```README.md``` located in the master branch.

In summary, the steps are:

1. Download rice reference assembly
2. Edit ```config.yaml``` so that:
  * Workflow recognizes rice assembly instead of potato
  * Workflow uses a rice bed file as a control for making dosage plots
3. Run ```scripts/units_from_SRA_RunTable.sh``` from this directory.
4. Run ```scripts/prefetch_rice.sh``` from this directory. Achtung, this uses UCD server modules to get the ascp executable!
5. Run ```scripts/parallel_fastq_dump.sh``` from this directory
6. Run workflow:

```
snakemake -s 1_init_genome_fofn.snakes --cores 8
snakemake -s 2_fastq_to_dosage_plot.snakes --cores 8
```

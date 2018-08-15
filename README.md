snakemake_chrom_dosage, a karyotyping tool

## Disclaimer

Kirk Amundson, 2018
This work is the property of the UC Davis Genome Center - Comai Lab

Use at your own risk.
We cannot guarantee support, but email kramundson at ucdavis dot edu
All information obtained/inferred with these scripts is without any implied warranty of fitness for any purpose or use whatsoever.

## Summary

This Snakemake workflow outputs read coverage by bin across a reference sequence, using
raw fastq reads in the subdirectory ```data/reads``` and a table of sample information,
```units.tsv``` as input. The bam file corresponding to control sample (including the path)
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

2. The last step in the workflow, plot generation, was previously done in JMP. For large
population datasets, reading these tables into memory became prohibitive. Reads are now
counted in non-overlapping bins of a user-specified size using bedtools. The bin size can
be changed by editing ```config.yaml```.

This workflow isn't nearly as user-friendly as bin-by-sam, but it is faster.

## Usage

Potato is used as an example.
s
1. Clone this repository

```
git clone https://github.com/kramundson/snakemake_chrom_dosage
cd snakemake_chrom_dosage
```

2. Install conda

Often, visitors don't have a home folder on the genome center clusters. In other words,
the folder ```home/yourUserName``` doesn't exist. If this is the case for you, skip to step
4 and see [here][./conda_alt_install.md] for a workaround.

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# review license
# accept license
# accept or change home location
# yes to placing it in your path

source $HOME/.bashrc
```

3. Build environment using included file ```environment.yaml```

```
conda env create --name dosage -f environment.yaml

# follow prompts if necessary
```

4. Symlink or copy your uninterleaved reads to be analyzed to the folder ```./data/reads```

5. Edit ```units.tsv``` to suit your needs.

units.tsv is a 5-column tab delimited text file. Each column specifies the following:

sample: unique biological sample
units: unique combination of biolgical sample, library and sequencing run
fq1: name of file corresponding to the sample-unit combination located in ```data/reads```
fq2: always NaN for this analysis (left in for paired-end data analyses, that takes the same type of table)
parhap: a portmenteau of parent-haploid, basically set to "mother" for control sample and "haploid" to test samples

6. Edit ```config.yaml``` to suit your needs.

The main thing to change is the name of the control file for the relative coverage analysis.
The leading path will always be ```data/bedtools_covereage/```, and the file name will be
the sample identifier (from units.tsv) for the control sample with a .bed extension

Command line parameters can also be changed here, if desired.

7. Optional: edit ```cluster.yaml``` to suit your needs (advanced users only).
The current specs should work for most applications.

8. Run the first part of the workflow.

You can (and should) do a Snakemake "dry run" to see that the workflow is set up correctly.

```
snakemake -s init_genome.snakes -np
```

To run it for real, remove the -np flag.

```
snakemake -s init_genome.snakes
```

This should take ~20 minutes for potato

9. Test the remaining steps with a dry run, then run the remaining steps of the workflow
on the cluster

```
snakemake -s Snakefile -np
sbatch runSnakes.slurm
```

Time from raw reads to dosage plots is ~1 day.

## Output

For each sample, one plot showing relative coverage values per bin, one table with the
raw read counts and relative coverage values for each bin. Currently, intermediate files
are kept (note, this takes up a lot of disk space).

Notes:

* By analysis-ready, I really mean that the chloroplast and mitochondrion assemblies
of the DM1-3 reference genotype are added to the standard DM1-3 v4.04 assembly, with some
necessary cleanup of the FASTA headers. I found that including the chloroplast during the
mapping, and then filtering for alignments of adequate mapping quality cleans up the plots
substantially.
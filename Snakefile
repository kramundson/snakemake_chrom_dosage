# Run init_genome.snakes before running this Snakefile
# This Snakefile handles dataset-specific analysis, assuming that reference genome
# file dependencies and file of output filenames have been made using init_genome.snakes

import re
import pandas as pd
shell.executable('bash')

configfile: "config.yaml"

units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

def is_single_end(sample,unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_fastq(wildcards):
    return "data/reads/"+units.loc[(wildcards.sample, wildcards.unit), ['fq1', 'fq2']].dropna()

def get_fastq_test(wildcards):
    tmp = units.loc[units['parhap'] != 'mother']
    return "data/reads/"+tmp.loc[(wildcards.sample,wildcards.unit), ['fq1','fq2']].dropna()

def get_fastq_control(wildcards):
    tmp = units.loc[units['parhap'] == 'mother']
    return "data/reads/"+tmp.loc[(wildcards.sample,wildcards.unit), ['fq1','fq2']].dropna()

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample should be aligned as such
        return expand("data/trimmed/{sample}-{unit}-{group}.fastq.gz",
            group=[1,2], **wildcards)
    # single end sample
    return "data/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

rule all:
    input: # TODO gather step for dosage plots, probably need gatekeeper pseudorules
        config["targets"]

include: "rules/dosage_plot_fofn.rules"
include: "rules/dosage_plots.rules"
# include: "rules/dosage_plot_population.rules"
include: "rules/bedtools_coverage.rules"
include: "rules/samtools_index.rules"
include: "rules/bam_mapqual_filter.rules"
include: "rules/samtools_merge.rules"
include: "rules/mark_duplicates.rules"
include: "rules/align.rules"
include: "rules/cutadapt_pe.rules"
include: "rules/cutadapt.rules"
include: "rules/cutadapt_hardtrim.rules"
include: "rules/get_SRA_reads.rules"

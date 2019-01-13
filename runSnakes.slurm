#!/bin/bash
#
#SBATCH --job-name=snake_mapping
#SBATCH -c 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=1G # Memory pool for all cores in MB (see also --mem-per-cpu)
#SBATCH --time 3-00:00:00 # 3 day runtime
#SBATCH -p production # Partition to submit to
#SBATCH -o pipeline_test.out # File to which STDOUT will be written
#SBATCH -e pipeline_test.err # File to which STDERR will be written

date
hostname
# cd /share/comailab/kramundson/snakemake_chrom_dosage # modify path as needed

source activate dosage

snakemake -j 7 -s 1_init_genome_fofn.snakes --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time} -c {cluster.c} --mem-per-cpu {cluster.mempercpu}" -k -w 120
snakemake -j 999 -s 2_fastq_to_dosage_plot.snakes --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time} -c {cluster.c} --mem-per-cpu {cluster.mempercpu}" -k -w 120

source deactivate

date

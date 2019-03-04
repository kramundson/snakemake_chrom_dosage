#!/bin/bash

# putting prefetch back in, but outside of snakemake workflow for now
# todo, how to implement this robustly when ascp executable and asperaweb_id_dsa.openssh could be anywhere?
# note, I altered the default prefetch path so reads would go right to this folder. Not ideal.
module load aspera-connect/3.5.1
cut -f 2 units.tsv | grep "^SRR" | xargs prefetch --ascp-path '/software/aspera-connect/3.5.1/static/bin/ascp|/software/aspera-connect/3.5.1/static/etc/asperaweb_id_dsa.openssh' --max-size 150G

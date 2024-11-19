import sys
import argparse
import pysam
import logging
import array
import numpy as np
from timeit import default_timer as timer
import math
import os
from itertools import islice

parser = argparse.ArgumentParser(
                    prog='consensus_vcf.py',
                    description='Takes a consensus genome and consensus-coverted alternate genomes and makes a vcf for identifying species.')
                    
parser.add_argument("-c", "--consensus",     required=True, help="The N-way consensus genome")
parser.add_argument("-g", "--genomes",       required=True, help="Genomes in consensus coordinates, separated by commas")
parser.add_argument("-l", "--labels",        required=True, help="Column headers to use in the VCF file for these genomes, separated by commas")
parser.add_argument("-o", "--output-prefix", required=True, help="The prefix for the output files")

args = parser.parse_args()

genome_files   = args.genomes.split(",")
labels         = args.labels.split(",")
output         = args.output_prefix
consensus_file = args.consensus

consensus_genome = pysam.Fastafile(consensus_file)
genomes = []
for genome in genome_files:
  genomes.append(pysam.Fastafile(genome))

vcfh = pysam.VariantHeader()
for label in labels: vcfh.add_sample(label)
for chrom in consensus_genome.references:
  vcfh.add_meta('contig', items=[('ID',chrom)])

# Add GT to FORMAT in our VCF header
vcfh.add_meta('FORMAT', items=[('ID',"GT"), ('Number',1), ('Type','String'),
    ('Description','Genotype')])

# Open a file, "example.vcf" for writing, with our "vcfh" header
vcf = pysam.VariantFile(output +".vcf", "w", header=vcfh)

for chrom in consensus_genome.references:
  consensus_chr = consensus_genome.fetch(chrom).upper()
  
  other_chrs = []
  for genome in genomes: other_chrs.append(genome.fetch(chrom).upper())
  
  l = len(consensus_chr)
  
  for i in islice(range(l), 1, l-1):
    
    # only stop on snps
    if consensus_chr[i] != 'N': continue
    if consensus_chr[i-1] == 'N': continue
    if consensus_chr[i+1] == 'N': continue
    
    bases = []
    for this_chr in other_chrs:
      bases.append(this_chr[i])
    
    # only include snps with exactly two values
    ubases = list(set(bases))
    if not len(ubases)==2: continue
    
    r = vcf.new_record(contig=chrom, start=i, stop=i+1,
            alleles=(ubases[0],ubases[1]))
    
    for j in range(len(bases)):
      if bases[j] == ubases[0]:
        r.samples[labels[j]]['GT'] = (0,0)
      else:
        r.samples[labels[j]]['GT'] = (1,1)
        
    vcf.write(r)

vcf.close()




import sys
import argparse
import pysam
import logging
import array
import numpy as np
from timeit import default_timer as timer
import math
import os

parser = argparse.ArgumentParser(
                    prog='multiple_consensus.py',
                    description='Takes consensus genomes (in the same coordinate system) builds a new consensus.')

parser.add_argument("-c", "--consensus",     required=True, help="A consensus genome. Use argument multiple times for an N-way consensus.", action='append')
parser.add_argument("-o", "--output-prefix", required=True, help="The prefix for the output files")

args = parser.parse_args()

genome_files = args.consensus
output       = args.output_prefix

if len(genome_files) < 2:
  raise Exception("Error: You must specify at least two genomes with option --consensus (-c)")

# hold the first genome out as the consensus genome
consensus_genome = pysam.Fastafile(genome_files[0])
del genome_files[0]

# store the other genomes as a list
genomes = []
for genome in genome_files:
  genomes.append(pysam.Fastafile(genome))
  
  
ofile = open(output + "_consensus.fa", "w")
for chrom in consensus_genome.references:
  start = timer()
  print("Processing chromosome " + chrom, file=sys.stderr)
  
  write = True
  
  for genome in genomes:
    if write == False: continue
    
    if not chrom in genome.references:
      print("Warning: chromosome " + chrom + " not found in all genomes. " +
                      "Were these all produced using make_consensus.py using the" +
                      " same reference?", file=sys.stderr)
      write = False
      continue
    
    consensus_chr = list(consensus_genome.fetch(chrom).upper())
    this_chr      = genome.fetch(chrom).upper()
    
    chr_length = len(consensus_chr)
    if not len(consensus_chr)==chr_length:
      raise Exception("Error: chromosome " + chrom + " is not the same length." +
                      "Were these all produced using make_consensus.py using the" +
                      " same reference?")
    
    for i in range(chr_length):
      if consensus_chr[i] == 'N': continue
      if this_chr[i] == 'N': consensus_chr[i] = 'N'
    
  consensus_chr = ''.join(consensus_chr)
  end = timer()
  print("Completed in", round(end-start,2), "seconds\n", file=sys.stderr)
  if write: ofile.write(">" + chrom + "\n" + consensus_chr + "\n")
  

ofile.close()
  

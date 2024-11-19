import sys
import argparse
import pysam
import logging
import array
import numpy as np
from timeit import default_timer as timer
from lib import liftover_functions as lift
import math
import os

parser = argparse.ArgumentParser(
                    prog='make_consensus.py',
                    description='Builds a 2-way consensus genome in the coordinate system of the reference genome')

parser.add_argument("-r", "--reference",     required=True, help="The fasta sequence of the reference genome")
parser.add_argument("-a", "--alternate",     required=True, help="The fasta sequence of the alternate genome")
parser.add_argument("-c", "--chain",         required=True, help="The reference to alternate chain file")
parser.add_argument('-l', "--length",        required=False, type=int, default=6, help='The number of bases around indels to mask')
parser.add_argument("-o", "--output-prefix", required=True, help="The prefix for the output files")

args = parser.parse_args()

reference   = args.reference
alternate   = args.alternate
chainfile   = args.chain
mask_length = args.length
output      = args.output_prefix

start = timer()
print("Reading chain file",file=sys.stderr)

REFFASTA = pysam.Fastafile(reference)
ALTFASTA = pysam.Fastafile(alternate)

maps, tSizes, qSizes = lift.read_chain_file(chainfile, REFFASTA.references, ALTFASTA.references)

end = timer()
print("Completed in", round(end-start,2), "seconds\n", file=sys.stderr)

# This function does all the work. It compares each interval base-by-base
# and replaces N with the nucleotide only if the sequences match in each
# reference. It starts and ends before the begging/end of intervals so that 
# any sequence around a indel stays N

def update_interval_consensus(interval, ref_chr, consensus_chr, reference_chr, alternate_chr):
  ref_start  = interval.begin
  ref_end    = interval.end
  
  alt_chr    = interval.data[0]
  alt_start  = interval.data[1]
  alt_end    = interval.data[2]
  alt_strand = interval.data[3]
  
  ref_seq = REFFASTA.fetch(ref_chr, ref_start, ref_end).upper()
  alt_seq = ALTFASTA.fetch(alt_chr, alt_start, alt_end).upper()
  
  if alt_strand == "-":
    alt_seq = lift.revcomp_DNA(alt_seq)
  
  l = len(ref_seq)
  
  ndiff = 0
  nsame = 0
  
  for i in range(mask_length,l-mask_length):
      this_base = ref_seq[i]
      alt_base  = alt_seq[i]
      if this_base == alt_base:
          nsame += 1
          consensus_chr[i+ref_start] = this_base
          reference_chr[i+ref_start] = this_base
          alternate_chr[i+ref_start] = this_base
      else:
          ndiff +=1
          reference_chr[i+ref_start] = this_base
          alternate_chr[i+ref_start] = alt_base
          
          
  tot = ndiff + nsame
  
  # throw an error is there is less than 50% identity in any interval. 
  # it probably means we did something wrong
  
  if tot > 100:
    if nsame/tot < 0.5:
      print("Warning: The interval " + 
        ref_chr + ":" + str(ref_start) + "-" + str(ref_end) +
        " had " + str(round(100*nsame/tot)) + 
        " identity! There may be a problem with your sequence or chain files", file=sys.stderr)
      print(ref_seq)
      print(alt_seq)
  return 1
  
# I want the highest scoring chain to be the main chain, so I run through chains
# in reverse order of score. To handle situations where multiple chains map
# to a single area, I overwrite the consensus with "N" before processing the chain.

consensus_genome = {}
reference_genome = {}
alternate_genome = {}

for this_chr, this_len in tSizes.items():
  if "alt" in this_chr: continue
  
  start = timer()
  print("Processing chromosome " + this_chr, file=sys.stderr)
  
  consensus_genome[this_chr] = list("N" * this_len)
  reference_genome[this_chr] = list("N" * this_len)
  alternate_genome[this_chr] = list("N" * this_len)

  this_chains = lift.get_chr_chains(maps, this_chr)
  chains = sorted(this_chains, key=lambda chain: chain.data['score'])
  
  for this_chain in chains:
    consensus_genome[this_chr][this_chain.begin:this_chain.end] = list("N" * (this_chain.end - this_chain.begin))
    reference_genome[this_chr][this_chain.begin:this_chain.end] = list("N" * (this_chain.end - this_chain.begin))
    alternate_genome[this_chr][this_chain.begin:this_chain.end] = list("N" * (this_chain.end - this_chain.begin))
    
    for interval in this_chain.data['mapTree'].items():
      val = update_interval_consensus(interval, 
        this_chr, 
        consensus_genome[this_chr], 
        reference_genome[this_chr], 
        alternate_genome[this_chr])
  
  end = timer()
  print("Completed in", round(end-start,2), "seconds\n", file=sys.stderr)

  start = timer()
  print("Validating this chromosome", file=sys.stderr)
  ref_ns     = 0
  consensus_ns = 0
  ref_seq = REFFASTA.fetch(this_chr, 0, this_len).upper()

  if len(ref_seq) != len(consensus_genome[this_chr]):
    print("Error: reference sequence and consensus sequence were not the same length for " + this_chr, file=sys.stderr)
    break
  
  for i in range(this_len):
    if ref_seq[i] == 'N':
      ref_ns += 1
    if consensus_genome[this_chr][i] == 'N':
      consensus_ns += 1
      continue
    if consensus_genome[this_chr][i] != ref_seq[i]:
      print("Error: reference sequence and consensus sequence not identical for " + 
      this_chr + " " + str(i), file=sys.stderr)
  
  print("Reference " + this_chr + " is " + str(round(100*ref_ns/this_len)) + "% N", file=sys.stderr)
  print("Consensus " + this_chr + " is " + str(round(100*consensus_ns/this_len)) + "% N", file=sys.stderr)
  
  end = timer()
  print("Completed in", round(end-start,2), "seconds\n", file=sys.stderr)
    
consensus_out = {}
reference_out = {}
alternate_out = {}
for this_chr, this_len in tSizes.items():
  if "alt" in this_chr: continue
  consensus_out[this_chr] = ''.join(consensus_genome[this_chr])
  reference_out[this_chr] = ''.join(reference_genome[this_chr])
  alternate_out[this_chr] = ''.join(alternate_genome[this_chr])
  

ofile = open(output + "_consensus.fa", "w")
for this_chr, this_seq in consensus_out.items():
    ofile.write(">" + this_chr + "\n" + this_seq + "\n")
ofile.close()

ofile = open(output + "_reference.fa", "w")
for this_chr, this_seq in reference_out.items():
    ofile.write(">" + this_chr + "\n" + this_seq + "\n")
ofile.close()

  
ofile = open(output + "_alternate.fa", "w")
for this_chr, this_seq in alternate_out.items():
    ofile.write(">" + this_chr + "\n" + this_seq + "\n")
ofile.close()

  


# ConsensusGenomeTools
Comparative genomics involves quanitiative comparisons of molecular features across species. These comparisons are confounded by alignment bias. One solution to this problem is to align all data to a single consensus genome. In the consensus genome 'N' is used to mask single nucleotide differences, as well as windows around insertions and deletions. This removes the majority of alignment bias, and also simplifies many downstream tasks because it places data from multiple species into a single reference coordinate system. As far as I am aware, this method was first introduced by [He, Bammann, and Han 2014](https://doi.org/10.1261/rna.043075.113).

This package includes 3 tools for creating and working with consensus genome files:
1. *make_consensus.py* builds a 2-way consensus genome from a reference genome, alternate genome, and a reference-to-alternate UCSC chain file.
2. *multiple_consensus.py* takes multiple 2-way consensus genomes that use the same reference and builds an N-way consensus genome.
3. *consensus_vcf.py* takes a consensus genome and produces a VCF file of all single-nucleotide fixed differences between genome refereneces.

In the examples for each tool I will provide scripts that could be used to generate a human/chimpanzee/gorilla consensus genome as well as a VCF with fixed differences between the three. 

## Installation
These tools are all python scripts. These will all run in the "pysam" conda environment that is specified by the yml in the env folder. 

## Notes
Fasta files will need to be indexed in order to provide random-access. This means that files cannot be zipped with gzip. If you do want your files zipped, please use bgzip to zip your files instead. You do not need to create the indexes yourself as pysam will create them if they don't exist.

Make sure that your chain files and your fasta files have matching contig names. Typically this will be UCSC-style in order to match the UCSC chain files. These tools will print a warning for any missing contigs and will continue processing only the contigs with matching names. 

Some contigs may drop out because there is no chain mapping them between species. This means your final consensus genome may be missing some contigs. For closely related species, these should genrally be small and unlocalized contigs that are unlikely to contain annotated genes. Still, if you are building a reference using the consensus genome and a human annotation you may need to remove these contigs from your annotation.

## make_consensus.py
This takes two fasta genomes: a reference, and an alternate genome. It lifts all sequences in the reference to the alternate genome using the best chain, then compares each nucleotide. It will mask any mismatches with an N. It will also mask a region around insertions or deletions with a user-specified length of Ns (6 by default). It provides three outputs:
1. The consensus genome (with single nucleotide differences between species as 'N')
2. The consensus genome with single nucleotide differences replaced by the reference sequence
3. The consensus genome with single nucleotide differences replaced by the alternate sequence

Outputs 2 and 3 greatly simplify the later construction of a VCF.

```
usage: make_consensus.py [-h] -r REFERENCE -a ALTERNATE -c CHAIN [-l LENGTH] -o OUTPUT_PREFIX

Builds a 2-way consensus genome in the coordinate system of the reference genome

options:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        The fasta sequence of the reference genome
  -a ALTERNATE, --alternate ALTERNATE
                        The fasta sequence of the alternate genome
  -c CHAIN, --chain CHAIN
                        The reference to alternate chain file
  -l LENGTH, --length LENGTH
                        The number of bases around indels to mask (default 6)
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        The prefix for the output files
```
### example
```
make_consensus.py -r hg38.fa -a panTro6.fa -c hg38ToPanTro6.over.chain.gz -o hg38_panTro6
make_consensus.py -r hg38.fa -a gorGor6.fa -c hg38ToGorGor6.over.chain.gz -o hg38_gorGor6
```

## multiple_consensus.py
This tool takes multiple 2-way consensus genomes, built by make_consensus.py, and builds an N-way consensus genome. 
```
usage: multiple_consensus.py [-h] -c CONSENSUS -o OUTPUT_PREFIX

Takes consensus genomes (in the same coordinate system) builds a new consensus.

options:
  -h, --help            show this help message and exit
  -c CONSENSUS, --consensus CONSENSUS
                        A consensus genome. Use argument multiple times for an N-way consensus.
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        The prefix for the output files
```
### example
```
multiple_consensus.py -c hg38_panTro6_consensus.fa -c hg38_gorGor6_consensus.fa -o hg38_panTro6_gorGor6
```

## consensus_vcf.py
This builds a vcf of fixed single-nucleotide differences between species. This is useful for identifiying species in multiplexed 'zoo' experiments in seingle cell data. After alignment to a consensus genome, the vcf can be used with tools like demuxlet or vireo to identify species. At least between human/chimpanzee/gorilla this results in a substantial number of SNPs, so it is best to filter the VCF to a reduced set. For instance taking only exonic SNPs.

```
usage: consensus_vcf.py [-h] -c CONSENSUS -g GENOMES -l LABELS -o OUTPUT_PREFIX

Takes a consensus genome and consensus-coverted alternate genomes and makes a vcf for identifying species.

options:
  -h, --help            show this help message and exit
  -c CONSENSUS, --consensus CONSENSUS
                        The N-way consensus genome
  -g GENOMES, --genomes GENOMES
                        Genomes in consensus coordinates, separated by commas
  -l LABELS, --labels LABELS
                        Column headers to use in the VCF file for these genomes, separated by commas
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        The prefix for the output files
```
### example
```
consensus_vcf.py -c hg38_panTro6_gorGor6_consensus.fa -g hg38_panTro6_reference.fa,hg38_panTro6_alternate.fa,hg38_gorGor6_alternate.fa -l human,chimp,gorilla -o human_chimp_gorilla
```

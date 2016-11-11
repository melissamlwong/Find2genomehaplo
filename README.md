
This program aims to identify two haplotype sequences within a genome assembly guided by a reference genome. A common problem in genome assembly of heterozygous genome using Illumina short reads is two haploid genome were assembled separately. The detection of these sequences is very important and interesting to some biologists, e.g. studying the origins, mutations and structural variations between the two haploid genomes especially from a hybrid species and investigating X-linked and Y-linked genes etc. 

Usage
usage: find2genomehaplo.py [-h] infile ref gff

Identify two haplotypes in a higly heterozygous genome assembly using
reference genome

positional arguments:
  infile      input fasta file
  ref         masked reference genome
  gff         gff3 file

optional arguments:
  -h, --help  show this help message and exit

##Detect two haplotypes in genome assembly using reference genome

This program aims to identify two haplotype sequences within a genome assembly guided by a reference genome. A common problem in genome assembly of heterozygous genome using Illumina short reads is two haploid genome were assembled separately. The detection of these sequences is very important and interesting to some biologists, e.g. studying the origins, mutations and structural variations between the two haploid genomes especially from a hybrid species and investigating X-linked and Y-linked genes etc. 

####Workflow:
1. Identify scaffolds containing exons from a reference genome using blastn
2. blast-2-sequences the top 2 best hit for each exon
3. Identify syntenic blocks using DAGchainer
4. Filter duplicated genes and genes with no exons

##Getting started

####Install the following dependencies:

* Executable for [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* Executable for [makeblastdb](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* Executable for [dagchainer](https://sourceforge.net/projects/dagchainer/)

*Note: Full path must be specified in the script if the executable is not found in the path

##Usage

```bash
usage: find2genomehaplo.py [-h] infile ref gff

Identify two haplotypes in a higly heterozygous genome assembly using
reference genome

positional arguments:
   infile      input fasta file
   ref         masked reference genome
   gff         gff3 file
optional arguments:
  -h, --help  show this help message and exit
```

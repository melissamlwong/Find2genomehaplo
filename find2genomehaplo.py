#!/usr/bin/env python

#This script is written by Melissa M.L. Wong (Melissa.Wong@unil.ch) to identify two differentiated haplotypes in draft genome assembly using exon-guided approach

import sys
import os
import commands
import re
import argparse
import subprocess
from os.path import splitext
from collections import defaultdict
from itertools import groupby, izip

##### Executables and parameters - Edit here! ##### 
dagchainer_fullpath = "/home/melissamlwong/softwares/DAGCHAINER/run_DAG_chainer.pl"
blastn_fullpath = "blastn"
makeblastdb_fullpath = "makeblastdb"
num_threads = 24 #argparse
filter_id = 90 #argparse
#Dagchainer option g, D, A
##################################################

def optparser():
    """Option parser"""
    p = argparse.ArgumentParser(description='Identify two haplotypes in a higly heterozygous genome assembly using reference genome')
    p.add_argument("infile", nargs='?', type=argparse.FileType('r'),default=sys.stdin,help="input fasta file")
    p.add_argument("ref", nargs='?', type=argparse.FileType('r'),default=sys.stdin,help="masked reference genome")
    p.add_argument("gff", nargs='?', type=argparse.FileType('r'),default=sys.stdin,help="gff3 file")
    args = p.parse_args()
    return args


def check_exefiles():
    print "*** Step 1 of 6 : Checking executables and input files ..."
    with open(sys.argv[-3],'r') as f:
        for line in f:
            if re.findall("|",line):
                sys.stderr.write("  * Error : Input fasta file must not contain \"|\" symbol!!!\n")
                Exit = True
    if len(sys.argv) == 4:
        print "  * Enough arguments in command line"
    else:
        sys.stderr.write("  * Error : Not enough arguments in command line!!!\n" )
        Exit = True
    if commands.getoutput("if [ ! -z $(command -v %s) ] && [ ! -z $(command -v %s) ] && [ ! -z $(command -v %s) ]; then echo True; else echo False; fi" % 
        (dagchainer_fullpath, blastn_fullpath, makeblastdb_fullpath)) == "True":
        print "  * The following executables are working:\n    1) %s\n    2) %s\n    3) %s" % (blastn_fullpath, makeblastdb_fullpath, dagchainer_fullpath)
    else:
        sys.stderr.write("  * Error : One or more executables do not exist!!!\n     Please insert the correct path in the script\n")
        Exit = True
    if os.path.isfile(sys.argv[-1]) and os.path.isfile(sys.argv[-2]) and os.path.isfile(sys.argv[-3]):
        print "  * All files exist ..."
    else:
        sys.stderr.write("  * Error : One or more files do not exist!!!\n" )
        Exit = True
    if Exit:
        sys.stderr.write("*** Exiting ***\n%s\n" % p.print_help())
        sys.exit()

##### Part 1 - Selection #####
def exonpos(gff):
    """Extract CDS/UTR positions from gff files and get their fasta sequence"""
    exon_dict = defaultdict(dict) #use a nested dictionary
    for line in gff:
        try:
            chrom,source,feature,start,stop,score,strand,frame,name = line.split("\t")
        except:
            continue
        else:
            if feature == "CDS" or feature == "five_prime_UTR" or feature == "three_prime_UTR":
                id_name = name.split(";")[0].split("=")[1]
                new = id_name.replace("five_prime_UTR","5UTR").replace("three_prime_UTR","3UTR")
                exon_dict[chrom][new]=[start,stop]
    if exon_dict:
        return exon_dict
    else:
        sys.stderr.write("  * Error in parsing gff3 file. Please check formatting.\n*** Exiting ***")
        sys.exit()

def extract_exon(ref,extract):
    """Extract exon sequences when given a nested dictionary"""
    fh = open(ref,'r')
    exon = open("exon.fa",'w')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        if header in extract.keys():
            for i in extract[header].keys():
                start,stop = extract[header][i]
                length = int(stop)-int(start)+1
                if length >=22: #The recommended minimum sequence length is 22 bases
                    exon.write(">%s\n%s\n" % (i,seq[int(start)-1:int(stop)]))
                    exon.flush()
                    #count +=1
    #print "Extracted a total of %s exonic sequence (min 22bp) from gff file" % count #count is missing

def extract_input(infile):
    """Extract input sequences when given a nested list"""
    unique_pair = []
    allseq_dict=defaultdict()
    with open("unique_pair.dat",'r') as u:
        for line in u:
            x,y=line.strip("\n").split("\t")
            unique_pair.append([x,y])
    l=[x for b in unique_pair for x in b]
    allseq = set(l)
    fh = open(infile,'r')
    out1 = open("out1.fa",'w')
    out2 = open("out2.fa",'w')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        if header in allseq:
            allseq_dict[header]=seq
    for h1,h2 in unique_pair:
            out1.write(">%s\t%s\n" % (h1,allseq_dict[h1]))
            out2.write(">%s\t%s\n" % (h2,allseq_dict[h2]))

def top2blastn(inputfile,num_threads,exe1,exe2):
    """blastn"""
    top2_dict = defaultdict(list)
    unique_pair = []
    pair=[]
    runblastn = [exe1,"-task","blastn","-db","ref","-query","exon_frogref.fa","-evalue","1E-5","-outfmt","6","-out","allexon.blast","-num_threads",num_threads,"-max_hsps","1"]
    print "  * Building blastn DB for input fasta file..."
    subprocess.call([exe2,"-in",inputfile,"-dbtype","nucl","-out","ref"])
    print "\n\n  * Blastn exon sequences to input fasta file...\n"
    subprocess.call(runblastn)
    print "\n\n  * Blastn completed\n"
    top2 = open("top2_blast.tagset",'w') 
    with open("allexon.blast") as f:
        for line in f:
            rseq,qseq,nucid,ali_l,ali_s,ali_gap,rstart,rstop,qstart,qstop,evalue,score = line.split("\t")
            if len(top2_dict[rseq])==0:
                top2_dict[rseq].append(qseq)
            elif len(top2_dict[rseq])==1:
                if qseq != "".join(top2_dict[rseq]):
                    top2.write("%s\t%s\t%s\n" % (rseq,"".join(top2_dict[rseq]),qseq))
                    top2_dict[rseq].append(qseq)
                    name = "\t".join(top2_dict[rseq])
                    pair.append(name)
                else:
                    top2_dict[rseq].append("Same")
    unique = open("unique_pair.dat",'w')
    for i in list(set(pair)):
        first,second = i.split("\t")
        reverse = "%s\t%s" % (second,first)
        if not i in unique_pair and not reverse in unique_pair:
            unique_pair.append(i)
            unique.write("%s\n" % i)

##### Part 2 - Detect synteny between any two related sequences#####
def blast2seq(filter_id,exe):
    """blast 2 sequences for each top2 pair and convert blast results to dagchainer input format"""
    f1 = open('out1.fa','r')
    f2 = open('out2.fa','r')
    blast = open(basename + "_b2seq.blast","w")
    dag = open(basename + "_dagchainer.dat","w")
    count = 0
    nohit=0
    for l1,l2 in izip(f1, f2):
        count+=1
        subject = open('subject.fa','w+')
        query = open('query.fa','w+')
        h1,s1=l1.split("\t")
        h2,s2=l2.split("\t")
    	query.write("%s\n%s\n" % (h1,s1))
    	subject.write("%s\n%s\n" % (h2,s2))
        subject.close()
        query.close()
    	p=subprocess.Popen([exe,"-task","blastn","-subject","subject.fa","-query","query.fa","-evalue","1E-5","-outfmt","6"],stdout=subprocess.PIPE)
        hits=0
        for line in p.stdout.readlines():
            if len(line.split("\t")) == 12:
                rseq,qseq,nucid,ali_l,ali_s,ali_gap,rstart,rstop,qstart,qstop,evalue,score = line.split("\t")
                blast.write("%s" % line)
                hits+=1
                if nucid >= filter_id:
                    dag.write("%s\t%s_1\t%s\t%s\t%s\t%s_2\t%s\t%s\t%s\n" % (rseq,count,rstart,rstop,qseq,count,qstart,qstop,evalue))
        if hits == 0: nohit+=1
    print "  * %s out of %s pairs (%.1f%) has no blastn hit" % (nohit,count,nohit/count)
    f1.close()
    f2.close()
    blast.close()
    dag.close()

def dagchainer(exe):
    infile = basename + "_dagchainer.dat"
    subprocess.call([exe,"-i",infile,"-s","-I","-g","1000","-D","10000","-A","3"])

if __name__ == '__main__':
    print "*** Initiating program ...\n"
    args = optparser()
    print "  * Status ok. Proceeding to Step 2\n"
    basename=splitext(sys.argv[0])[0] #fullpath to script
    exon_dict = exonpos(args.gff)
    args.gff.close() 
    extract_exon(ref,exon_dict)
    print "*** Step 2 of 6 : Extracting exon sequences completed\n"
    top2blastn(inputfile,str(num_threads),blastn_fullpath,makeblastdb_fullpath)
    print "*** Step 3 of 6 : Blastn exon sequences to input genome completed\n"
    extract_input(inputfile)
    print "*** Step 4 of 6 : Extracting input sequences completed\n"
    blast2seq(filter_id,blastn_fullpath)
    print "*** Step 5 of 6 : Blast 2 sequences completed\n"
    dagchainer(dagchainer_fullpath)
    print "*** Step 6 of 6 : Dagchainer completed\n"
    print "*** Job completed"
    args.ref.close() 
    args.infile.close()
    os.remove("subject.fa")
    os.remove("query.fa")
    os.remove("out1.fa")
    os.remove("out2.fa")

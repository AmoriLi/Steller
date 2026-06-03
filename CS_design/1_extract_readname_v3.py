import sys
import os
import argparse
import glob
import itertools
import gzip
from Bio import SeqIO
from scipy.spatial.distance import hamming
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import PairwiseAligner
#import pandas as pd
import time
from Bio.SeqRecord import SeqRecord

wd = sys.argv[1] #"NPC_project/BMK/cs_design/240122_batch2/PT-CS3/amplicon/highBC"
print(wd)
input = glob.glob(os.path.join(wd,"fastp")+"/"+"*1.fq.gz")[0]
print(input)
outdir = os.path.join(wd,"split_fastq")
ft1 = os.path.join(outdir,"cs3")
ft2 = os.path.join(outdir,"polyT")
ft3 = os.path.join(outdir,"other")
if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(ft1):
    os.makedirs(ft1)
if not os.path.exists(ft2):
    os.makedirs(ft2) 
if not os.path.exists(ft3):
    os.makedirs(ft3) 
os.chdir(wd)
os.getcwd()

cs3_probe = Seq("CTCAATAAAGCTTGCCTTGA")
cs3_probe_rev = cs3_probe.reverse_complement()
cs3_normal_rev = Seq("GCACTCAAGGCAAGCTTTATTGAG")
polyt = Seq("TTTTTTTTTTTTTTT")

# define alignment penalty
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = -0.5
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.5
aligner.target_end_gap_score = -0.5
aligner.query_end_gap_score = 0.0

### read fastq.gz
with gzip.open(input,"rt") as handle:
    seq_name = []
    seq = []
    for seq_record in SeqIO.parse(handle,"fastq"):
        seq_name.append(seq_record.id)
        seq.append(seq_record.seq)

polyt_read_name = []
polyt_read = []
polyt_cs3_read = []
polyt_cs3_read_name = []
#polyt_endindex = []
cs3_read_name = []
cs3_read = []
#cs3_endindex = []
other_read_name = []
other_read = []
    
    
for (read_name,read) in zip(seq_name,seq):
    # first estimate if polydT existed
    alg1 = aligner.align(read,polyt)
    if len(alg1) > 0:
        alg1 = aligner.align(read,polyt)[-1]
        alg1_score = alg1.score
        tmp=alg1.indices[1]
        tmp=tmp.tolist()
        polyt_endindex=tmp.index(max(tmp))
        if alg1_score >= 13 and polyt_endindex <= 130:
            subread = read[polyt_endindex:polyt_endindex+23]
            if len(subread) > 0:
                alg2 = aligner.align(subread,cs3_probe_rev)
                if len(alg2)>0:
                    alg2 = aligner.align(subread,cs3_probe_rev)[0]
                    alg2_score = alg2.score
                    if alg2_score >= 14:
                            #rev_read_name = read_name + " " + s
                        cs3_read_name.append(read_name)
                        cs3_read.append(read)
                    else:
                        polyt_read_name.append(read_name)
                        polyt_read.append(read)
                else:
                    polyt_read_name.append(read_name)
                    polyt_read.append(read)
            else:
            #m += 1 #read end with
                polyt_read_name.append(read_name)
                polyt_read.append(read)
        else:
            other_read_name.append(read_name)
            other_read.append(read)
    else:
        other_read_name.append(read_name)
        other_read.append(read)

with open(os.path.join(ft1,"cs3_readname.txt"), "w") as f:    
    for i in cs3_read_name:
        f.write(i + "\n")

with open(os.path.join(ft2,"polyT_readname.txt"), "w") as f:
    for i in polyt_read_name:
        f.write(i + "\n")

with open(os.path.join(ft3,"other_readname.txt"), "w") as f:
    for i in other_read_name:
        f.write(i + "\n")

print("cs3 read number is: " + str(len(cs3_read_name)))
print("polyT read number is: " + str(len(polyt_read_name)))
print("other read number is: " + str(len(other_read_name)))

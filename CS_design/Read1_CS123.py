## bulk amplicon of 4T1 tumor total RNA with high CloneBC transfection, which was captured by read1-cs1/2/3 bead oligos.
import sys
import os
import argparse
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
import glob

indir = sys.argv[1] #"~/NPC_project/BMK/cs_design/240122_batch2/read1-cs123/amplicon/highBC"
input = os.path.join(indir,"fastq/Unknown_YF-1-240112-BD-0115_good_1.fq.gz")
outdir = os.path.join(indir,"split_fastq")
ft1 = os.path.join(outdir,"cs1")
ft2 = os.path.join(outdir,"cs2")
ft3 = os.path.join(outdir,"cs3")
ft4 = os.path.join(outdir,"other")
if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(ft1):
    os.makedirs(ft1)
if not os.path.exists(ft2):
    os.makedirs(ft2) 
if not os.path.exists(ft3):
    os.makedirs(ft3) 
if not os.path.exists(ft4):
    os.makedirs(ft4)
os.chdir(indir)
os.getcwd()

cs1_probe_rev = Seq("GTCTTAAAGGTACTCTAG")
cs2_probe_rev = Seq("ATCTACAGCTGCCTTGTAAGTC")
cs3_probe_rev = Seq("TCAAGGCAAGCTTTATTGAG")

# define alignment penalty
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 1
aligner.mismatch_score = -0.5
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.5
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0

### read fastq.gz
with gzip.open(input,"rt") as handle:
    seq_name = []
    seq = []
    for seq_record in SeqIO.parse(handle,"fastq"):
        seq_name.append(seq_record.id)
        seq.append(seq_record.seq)

### align to polyt, cs3
cs1_read_name = []
cs1_read = []
cs2_read = []
cs2_read_name = []
#polyt_endindex = []
cs3_read_name = []
cs3_read = []
#cs3_endindex = []
other_read_name = []
other_read = []
for (read_name,read) in zip(seq_name,seq):
    # first estimate if polydT existed
    subread = read[0:25]
    alg1 = aligner.align(subread,cs1_probe_rev)
    if len(alg1) > 0:
        alg1 = aligner.align(subread,cs1_probe_rev)[0]
        alg1_score = alg1.score
        if alg1_score >= 16.5:
            cs1_read_name.append(read_name)
            cs1_read.append(read)
        else:
            alg2 = aligner.align(subread,cs2_probe_rev)
            if len(alg2)>0:
                alg2 = aligner.align(subread,cs2_probe_rev)[0]
                alg2_score = alg2.score
                if alg2_score >= 20.5:
                   cs2_read_name.append(read_name)
                   cs2_read.append(read)    
                else:
                    alg3 = aligner.align(subread,cs3_probe_rev)
                    if len(alg3) > 0:
                        alg3 = aligner.align(subread,cs3_probe_rev)[0]
                        alg3_score = alg3.score
                        if alg3_score >= 18.5:
                            #rev_read_name = read_name + " " + s
                            cs3_read_name.append(read_name)
                            cs3_read.append(read)
                            #cs3_endindex.append(alg1_index)
                        else:
                            other_read_name.append(read_name)
                            other_read.append(read)
                    else:
                        other_read_name.append(read_name)
                        other_read.append(read)
            else:
                other_read_name.append(read_name)
                other_read.append(read)
    else:
        other_read_name.append(read_name)
        other_read.append(read)

with open(os.path.join(ft1,"cs1_readname.txt"), "w") as f:
    for i in cs1_read_name:
        f.write(i + "\n")

with open(os.path.join(ft2,"cs2_readname.txt"), "w") as f:
    for i in cs2_read_name:
        f.write(i + "\n")

with open(os.path.join(ft3,"cs3_readname.txt"), "w") as f:
    for i in cs3_read_name:
        f.write(i + "\n")

with open(os.path.join(ft4,"other_readname.txt"), "w") as f:
    for i in other_read_name:
        f.write(i + "\n")

##2026-03
#AA, BB, AB, BA sampling PCR
import os
import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import re

def classify_reads_logic(input_gz, pt_seq, pt_rc,f_a, r_a, f_b, r_b):
    #
    files = {
        "AA": open("AA_read.txt", "w"),
        "BB": open("BB_read.txt", "w"),
        "AB": open("AB_read.txt", "w"),
        "BA": open("BA_read.txt", "w"),
        "None": open("trucated_read.txt", "w"),
        "Discarded": open("unmapped_read.txt", "w")
    }
    stats = {"AA": 0, "BB": 0,"AB": 0, "BA": 0, "None": 0, "Discarded": 0}
    #
    re_fwd = re.compile(pt_seq)
    re_rev = re.compile(pt_rc)

    print(f"Processing: {input_gz} ...")

    with open(input_gz, "r") as f:
        while True:
            #
            line1 = f.readline() # ID
            if not line1: break  # 
            line2 = f.readline().strip().upper() # 
            line3 = f.readline() # +
            line4 = f.readline() #
            seq_str=line2
            final_seq=None
            target_sub_seq = None
            # 
            match_fwd = re_fwd.search(seq_str)
            if match_fwd:
                final_seq = seq_str
                target_sub_seq = match_fwd.group()
            else:
                # 
                match_rev = re_rev.search(seq_str)
                if match_rev:
                    # 
                    final_seq = str(Seq(seq_str).reverse_complement()).upper()
                    # 
                    match_sub = re_fwd.search(final_seq)
                    if match_sub:
                        target_sub_seq = match_sub.group()

            # ------
            if final_seq and target_sub_seq:
                #
                if f_a in final_seq and r_a in final_seq:
                    files["AA"].write(target_sub_seq + "\n")
                    stats["AA"] += 1
                #
                elif f_b.upper() in final_seq and r_b.upper() in final_seq:
                    files["BB"].write(target_sub_seq + "\n")
                    stats["BB"] += 1
                #
                elif f_a.upper() in final_seq and r_b.upper() in final_seq:
                    files["AB"].write(target_sub_seq + "\n")
                    stats["AB"] += 1  
                elif f_b.upper() in final_seq and r_a.upper() in final_seq:
                    files["BA"].write(target_sub_seq + "\n")
                    stats["BA"] += 1 
                #
                else:
                    files["None"].write(target_sub_seq + "\n")
                    stats["None"] += 1
            else:
                files["Discarded"].write(seq_str + "\n")
                stats["Discarded"] += 1

    #
    for f in files.values():
        f.close()

    #
    print("\n--- Done ---")
    print(f"  └─  AA read: {stats['AA']}")
    print(f"  └─ BB read: {stats['BB']}")
    print(f"  └─  AB read: {stats['AB']}")
    print(f"  └─ BA read: {stats['BA']}")
    print(f"  └─ trucated read: {stats['None']}")
    print(f"unmmaped read: {stats['Discarded']}")
    
wd=sys.argv[1]
input=sys.argv[2]
pattern=r"[ATCG]{4}CT[ATCG]{4}AC[ATCG]{4}TC[ATCG]{4}GT[ATCG]{4}TG[ATCG]{4}CA[ATCG]{4}"
pattern_rv=r"[ATCG]{4}TG[ATCG]{4}CA[ATCG]{4}AC[ATCG]{4}GA[ATCG]{4}GT[ATCG]{4}AG[ATCG]{4}"
#
try:
    regex_pt = re.compile(pattern, re.IGNORECASE)
    print("successful！")
except re.error as e:
    print(f"Failed，position is {e.pos}: {e.msg}")
#pattern_rv="[ATCG]{4}TG[ATCG]{4}CA[ATCG]{4}AC[ATCG]{4}GA[ATCG]{4}GT[ATCG]{4}AG[ATCG]{4}"

af="GGCCACGCT" #GGCCACGCTCAGACGAGTCGGATCTCCCT
ar="AGCGCAGCC" #GACCAATGACTTACAAGGCAGCAGCGCAGCC
bf="GGCCGTGAG" #GGCCGTGAGCAGACGAGTCGGATCTCCCT
br="GAACTTGCC" #GACCAATGACTTACAAGGCAGCGAACTTGCC

# ------
classify_reads_logic(
    input_gz = input, #PE merged and fastqc filtered reads
    pt_seq   = pattern,      # barcode PT sequence
    pt_rc = pattern_rv,
    f_a      = af,       # primer A forward 
    r_a      = ar,       # primer A reverse-complement
    f_b      = bf,       # primer B forward
    r_b      = br       # primer B reverse-complement
)
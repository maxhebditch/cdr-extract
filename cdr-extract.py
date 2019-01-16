#!/usr/bin/env python
import re
import argparse

results = {}

class antibody:
    def __init__(self, name):
        self.name = name
        self.sequence = sequence

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

fasta_file = args.filename
chain = "Heavy"

with open(fasta_file,"r") as fasta:
    ab_name = ""
    sequence = ""
    sequence_found = False
    for line in fasta:
        if line.startswith(">"):
            ab_name = line.replace("\n","")
            ab_name = ab_name.replace(">","")
        else:
            sequence = line.replace("\n","")
            sequence_found = True
        if ab_name != "" and sequence_found != False:
            results[ab_name] = antibody(ab_name)
            ab_name = ""
            sequence_found = False

def build_CDR_and_non_array(input_sequence,input_index):
    CDR_all_index = []
    for CDR in input_index:
        for i in range(CDR[0],CDR[1]):
            CDR_all_index.append(i)

    CDR_combined_sequence = []
    non_CDR_combined_sequence = []
    for i, aa in enumerate(input_sequence):
        if i in CDR_all_index:
            CDR_combined_sequence.append(aa)
        else:
            non_CDR_combined_sequence.append(aa)
    
    return CDR_combined_sequence, non_CDR_combined_sequence

def find_sub_list(sl,l):
    results=[]
    sll=len(sl)
    for ind in (i for i,e in enumerate(l) if e==sl[0]):
        if l[ind:ind+sll]==sl:
            results.append((ind,ind+sll-1))
    return results

def getindex(cdr,chain):
    indx = find_sub_list(cdr,chain)[0]
    start = indx[0]
    end = indx[1]+1
    return start, end

def CDRH_finder(sequence):
    CDRH_index = []
    aseq = list(sequence)
    #CDRH1
    skip=19
    tseq = aseq[skip:]
    first_C = tseq.index("C")
    four_after_first_C = skip+first_C+4
    start_CDRH1 = four_after_first_C
    rough_CDRH1 = aseq[start_CDRH1:start_CDRH1+13]
    term_W = rough_CDRH1[::-1].index("W")
    CDRH1_a = rough_CDRH1[:-term_W-1]
    pos = find_sub_list(CDRH1_a,aseq)
    if len(pos) > 1:
        print("multiple CDRH1 in sequence")
    else:
        pos = pos[0]
    CDRH1_start = pos[0]
    CDRH1_stop = pos[1]-1
    CDRH1 = "".join(CDRH1_a)
    CDRH1_start, CDRH1_end = getindex(CDRH1_a,aseq)
    if len(CDRH1) > 12 or len(CDRH1) < 10:
        print("CDRH1 not the right size")
        raise Exception("CDR1 not identified")
        
    #CDRH2
    CDRH2_start = CDRH1_stop + 16
    rough_CDRH2 = aseq[CDRH2_start:CDRH2_start+15]
    end_CDRH2 = rough_CDRH2[-5:]

    for aa in range(0,6):
        close_aa = aseq[CDRH2_start+10+aa]
        far_aa = aseq[CDRH2_start+10+aa+33]
        if far_aa == "C":
            CDRH2_stop = CDRH2_start+10+aa+4
            CDRH3_start = CDRH2_start+10+aa+33+3
    try:
        CDRH2_stop
    except:
        raise Exception("CDR2 not identified")

    CDRH2_a = aseq[CDRH2_start:CDRH2_stop]
    CDRH2 = "".join(CDRH2_a)
    CDRH2_start, CDRH2_end = getindex(CDRH2_a,aseq)
        
    #CDRH3
    rough_CDRH3 = aseq[CDRH3_start:]
    CDRH3_end_found = False
    multi_pos = find_sub_list(['W','G'],rough_CDRH3)
    if len(multi_pos) > 0:
        for pot_pos in multi_pos:
            if rough_CDRH3[pot_pos[1]+2] == "G":
                CDRH3_end_found = True
                end_CDRH3 = pot_pos[0]
    if CDRH3_end_found == False:
        for ind, aa in enumerate(rough_CDRH3):
            if aa == "W":
                if rough_CDRH3[ind+3] == "G":
                    CDRH3_end_found = True
                    end_CDRH3 = ind
   
    try:
        end_CDRH3
    except:
        raise Exception("CDR3 not identified")

    CDRH3_a = rough_CDRH3[0:end_CDRH3]
    CDRH3 = "".join(CDRH3_a)
    CDRH3_start, CDRH3_end = getindex(CDRH3_a,aseq)
    
    CDRH_index.append([CDRH1_start,CDRH1_end])
    CDRH_index.append([CDRH2_start,CDRH2_end])
    CDRH_index.append([CDRH3_start,CDRH3_end])

    return CDRH1, CDRH2, CDRH3, CDRH_index

working_list = []
for antibody in results:
    sequence = (results[antibody].sequence)
    if chain == "Heavy":
        try:
            CDRH1, CDRH2, CDRH3, CDRH_index = CDRH_finder(sequence)
        except:
            print("not all CDRs identified for %s" % antibody)
        else:
            working_list.append(antibody)
            VH_CDR_combined, VH_non_CDR_combined = build_CDR_and_non_array(sequence,CDRH_index)
            results[antibody].CDR1 = CDRH1
            results[antibody].CDR2 = CDRH2
            results[antibody].CDR3 = CDRH3
            results[antibody].CDR_combined = "".join(VH_CDR_combined)
            results[antibody].non_CDR_combined = "".join(VH_non_CDR_combined)
    
with open("all_extracted.fasta","w") as outFile:
    for antibody in working_list:
        outFile.write(">"+antibody+"_whole_sequence\n")
        outFile.write(results[antibody].sequence+"\n")
        outFile.write(">"+antibody+"_CDR1\n")
        outFile.write(results[antibody].CDR1+"\n")
        outFile.write(">"+antibody+"_CDR2\n")
        outFile.write(results[antibody].CDR2+"\n")
        outFile.write(">"+antibody+"_CDR3\n")
        outFile.write(results[antibody].CDR3+"\n")
        outFile.write(">"+antibody+"_non_CDR_combined\n")
        outFile.write(results[antibody].non_CDR_combined+"\n")
        outFile.write(">"+antibody+"_CDR_combined\n")
        outFile.write(results[antibody].CDR_combined+"\n")

#!/usr/bin/env python

import argparse
import matplotlib
import matplotlib.pyplot as plt


def get_arguments():
    parser = argparse.ArgumentParser(description="This script reads through a \
fastq file and generates per base quality scores and average quality score graphs.")
    parser.add_argument("-i", "--infile1", help="Takes a .fa file containing \
the first set of paired reads in",\
                        required=True, type=argparse.FileType('r'))
    parser.add_argument("-j", "--infile2", help="Takes a .fq file containing \
the first index set in",\
                        required=True, type=argparse.FileType('r'))
    parser.add_argument("-k", "--infile3", help="Takes a file containing \
the indexes for our samples",\
                        required=True, type=argparse.FileType('r'))
    parser.add_argument("-o", "--outfile", help="Name of the index table output file",\
                        required=True, type=argparse.FileType('w'))
    parser.add_argument("-p", "--outfile2", help="Name of the swapped distribution output file",\
                        required=True, type=argparse.FileType('w'))    
    return parser.parse_args()

args=get_arguments()

# Reverse complement function
def rev_comp(dna):
    complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
    return ''.join([complement[base] for base in dna[::-1]])

# List of our indexes
indexlist = []
# Dictionary of acceptable index pairs
indexpairs = {}
# Dictionary of index sequence and ID
indexreference={}

# Read through index file and store our indexes
with args.infile3 as indexes:
    # Get rid of header line
    indexes.readline()
    for line in indexes:
        # Split line by tabs
        splitline = line.strip().split("\t")
        # Add index to indexes list
        indexlist.append(splitline[4])
        # Make dictionary entry where key is "index index"
        indexpairs[(splitline[4] + " " + splitline[4])] = \
                                 indexpairs.get((splitline[4] + " " + splitline[4]), 0)
        # Add index sequence and index ID to dictionary
        indexreference[splitline[4]] = splitline[3]

# Also add index swapped, undetermined, and sequencing error
indexpairs["index swapped"] = indexpairs.get("index swapped", 0)
indexpairs["undetermined"] = indexpairs.get("undetermined", 0)
indexpairs["sequencing error"] = indexpairs.get("sequencing error", 0)

# Create dictionary to store all properly indexed and index swapped samples
indexcounts={}

# Read through file and add quality scores to 
with args.infile1 as fh1:
    with args.infile2 as fh2:
        for line in fh1:
            # Read 4 lines from R2
            headr2=line.strip()
            seqr2=fh1.readline().strip()
            plusr2=fh1.readline().strip()
            qualr2=fh1.readline().strip()

            # Read 4 lines from R3
            headr3=fh2.readline().strip()
            seqr3=fh2.readline().strip()
            plusr3=fh2.readline().strip()
            qualr3=fh2.readline().strip()

            # If either indexes contain N, add 1 to undetermined
            if ("N" in seqr2) or ("N" in seqr3):
                indexpairs["undetermined"]+=1
                    
            # If both one of the indexes is not in the index list
            # add to sequencing error
            elif (seqr2 not in indexlist) or (rev_comp(seqr3) not in indexlist):
                indexpairs["sequencing error"]+=1

            # If indexes aren't undetermined or sequence errors
            # Check for index swapping
            elif (seqr2 + " " + rev_comp(seqr3)) not in indexpairs:
                indexpairs["index swapped"]+=1
                # Also add index pairs to index counts dict
                i1 = indexreference[seqr2]
                i2 = indexreference[rev_comp(seqr3)]
                indexcombo = i1 + " " + i2
                indexcounts[indexcombo] = indexcounts.get(indexcombo, 0) + 1

            # Otherwise, they are real read pairs
            else:
                indexpairs[seqr2 + " " + rev_comp(seqr3)]+=1
                # Also add index pairs to index counts dict
                i1 = indexreference[seqr2]
                i2 = indexreference[rev_comp(seqr3)]
                indexcombo = i1 + " " + i2
                indexcounts[indexcombo] = indexcounts.get(indexcombo, 0) + 1

# Print out index statistics out as table
with args.outfile as outfile:
    print("Index Pair", "Count", "Percentage", sep="\t", file=outfile)
    for key in indexpairs:
        print(key, indexpairs[key], (indexpairs[key]/sum(indexpairs.values())*100), sep="\t", file=outfile)

# Print index distribution out as table
with args.outfile2 as outfile2:
    print("Index 1", "Index 2", "Count", sep="\t", file=outfile2)
    for key in indexcounts:
        print(key.split(" ")[0], key.split(" ")[1], indexcounts[key], \
              sep="\t", file=outfile2)

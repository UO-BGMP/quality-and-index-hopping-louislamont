#!/usr/bin/env python

import argparse
import matplotlib
import matplotlib.pyplot as plt


def get_arguments():
    parser = argparse.ArgumentParser(description="This script reads through a \
fastq file and generates per base quality scores and average quality score graphs.")
    parser.add_argument("-i", "--infile", help="Takes a .fa file in",\
                        required=True, type=argparse.FileType('r'))
    parser.add_argument("-o", "--outfile", help="Name of the output file",\
                        required=True, type=str)
    parser.add_argument("-p", "--outfile2", help="Name of the output file",\
                        required=True, type=str)
    parser.add_argument("-r", "--rname", help="Name/Number of the read",\
                        required=True, type=str)
    return parser.parse_args()

args=get_arguments()

def convert_phred(base):
    '''Converts ASCII character to Phred score. Phred score = ASCII-33'''
    return ord(base)-33

# Read through file and add quality scores to array by bp
# Also append average quality score for the line to 
with args.infile as fh:
    # Initialize line count and array to hold all average read means
    NR = 0
    allmeans=[]

    # Read through file
    for line in fh:
        # Set mean qual score of line at 0 and increment line counter
        meanline=0
        NR += 1
        line=line.strip("\n")
        # If we're on the first sequence line, create a list the length of that
        # line to hold our running total of mean scores
        if NR == 2:
            mean_scores=[0.0]*len(line)

        # If we're on a quality line
        if NR % 4 == 0:
            # Iterate through the each base of that line
            for j in range(len(line)):
                # Convert the phred score into a quality score,
                # then add that to the same position in the mean_scores list
                # also add that score to the mean line score
                mean_scores[j]+=convert_phred(line[j])
                meanline+=convert_phred(line[j])

            # When we've read the whole line, divide meanline by the length
            # to get the average quality score for the line and add that to
            # allmeans list
            meanline=meanline/len(line)
            allmeans.append(meanline)
        if NR % 1000000 == 0:
            print("Processing line: ", NR, sep="")
            
# After reading all the lines, divide by the number of quality lines (NR/4)
for i in range(len(mean_scores)):
    mean_scores[i]=mean_scores[i]/(NR/4)

# Make phred score/bp plot and output to file
plt.title("Avg Phred Score per base position, " + args.rname)
plt.xlabel("Base position", size="large")
plt.ylabel("Phred Score", size="large")
plt.ylim(15,42)
plt.plot(mean_scores)
plt.savefig(args.outfile)
plt.clf()

# Make quality score histogram and output to file
# use log scale on Y-axis
plt.title("Histogram of read quality, " + args.rname)
plt.xlabel("Average quality score across read", size="large")
plt.ylabel("log(# of reads)", size="large")
plt.hist(allmeans)
plt.yscale("log", nonposy="clip")
plt.savefig(args.outfile2)
plt.clf()

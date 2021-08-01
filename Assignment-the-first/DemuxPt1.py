import Bioinfo
import gzip
import argparse


def get_args():
    parser = argparse.ArgumentParser(
        description="Plots median quality scores for a given fastq file"
    )
    parser.add_argument("-l", "--length", help="filename", required=True, type=int)
    return parser.parse_args()

args = get_args()

def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    for i in range(args.length):
        lst.append(value)
    return lst
my_list: list = []
my_list = init_list(my_list)

file = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
def populate_list(file):
    """Takes fastq file and returns list with total quality score by index and total number of lines
    in the file as an int"""
    with gzip.open(file, "rt") as fh:
        ls = []
        init_list(ls)
        counter = 0
        for line in fh:
            line = line.strip()
            counter += 1
            if counter % 4 == 0:
                for index, char in enumerate(line):
                    ls[index] += Bioinfo.convert_phred(char)
        return ls, counter
my_list, num_lines = populate_list(file)

print("# Base Pair\tMean Quality Score")
for index, total in enumerate(my_list):
    my_list[index] = total/(num_lines/4)
    print(str(index) + "\t" + str(my_list[index]))

import matplotlib.pyplot as plt
plt.bar(x=range(args.length), height=my_list)
plt.ylabel("Average Quality Score")
plt.xlabel("Index in Seqeunce")
plt.savefig("Histogram_r4.png")
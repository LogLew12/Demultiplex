import argparse
import Bioinfo
import gzip
def get_args():
    parser = argparse.ArgumentParser(
        description="Takes 4 fastq files (2 read files and 2 index files) and one file containing known barcodes. Populates files with reads by index and by reads. Places swapped indexes and unknonwn indexes into into separate files as well."
    )
    parser.add_argument("-r1", "--read1", help="filename of the first read file", required=True)
    parser.add_argument("-r2", "--read2", help="filename of the second read file", required=True)
    parser.add_argument("-i1", "--index1", help="filename of the first index file", required=True)
    parser.add_argument("-i2", "--index2", help="filename of the second read file (assumed to be reverse complemented)", required=True)
    parser.add_argument("-b", "--barcodes", help="filename of file containing all expected barcodes", required=True)
    parser.add_argument("-c", "--cutoff", type=int, help="quaility score cutoff for index reads. Any indexes with an average q-score below this will be put in the unknown file.", required=True)
    return parser.parse_args()

args = get_args()

fh_r1 = gzip.open(args.read1, "rt")
fh_r2 = gzip.open(args.read2, "rt")
fh_i1 = gzip.open(args.index1, "rt")
fh_i2 = gzip.open(args.index2, "rt")
fh_bars = open(args.barcodes, "r")
# fh_r1 = open(args.read1, "r")
# fh_r2 = open(args.read2, "r")
# fh_i1 = open(args.index1, "r")
# fh_i2 = open(args.index2, "r")

#creates a dictionary of barcodes, with sequence as the key and its "name" (ex. B1) as value
#also opens two correspoding files for writing (one for read 1 and read 2)
expected_barcode = {}
files1 = {}
files2 = {}
for line in fh_bars:
    line = line.strip()
    name = line.split("\t")[3]
    barcode = line.split("\t")[4]
    if Bioinfo.validate_base_seq(barcode) == True:
        expected_barcode[barcode] = name
        f1 = open(f"{name}_R1_out.fastq", "w") 
        f2 = open(f"{name}_R2_out.fastq", "w") 
        files1[barcode] = f1
        files2[barcode] = f2

#openning swapped and unknown files
swap_r1 = open("swapped_R1_out.fastq", "w")
swap_r2 = open("swapped_R2_out.fastq", "w")
unknown_r1 = open("unknown_R1_out.fastq", "w")
unknown_r2 = open("unknown_R2_out.fastq", "w")
#I maybe could have used file.writelines() instead of this
def write_file(ls, out):
    '''takes a list and an output file and writes every item
    as a line in the file'''
    for i in ls:
        out.write(i + "\n")
    return None

def revcomp(dna):
    '''takes a string of DNA and returns the reverse complement'''
    dna = dna.upper()
    dna_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    table = dna.maketrans(dna_dict)
    new_dna = dna.translate(table)
    new_dna = new_dna[::-1]
    return new_dna

def check_qscore(qscore_line, cutoff):
    '''Checks the average score in a quaility score line. Takes a quality score line as a string, and a cutoff as an int.
     Returns True if the average is less than the given cutoff; Returns False if not. '''
    total = 0
    for char in qscore_line:
        total += Bioinfo.convert_phred(char)
    if total/len(qscore_line) < cutoff:
        return True
    return False
def up_counter(count_dict, key):
    '''takes a dictionary of counters. If the given key exists, it is increased by 1, otherwise it is initiated to 1'''
    try:
        count_dict[key] += 1
    except:
        count_dict[key] = 1
#initializing lists to hold each record
record_read1 = []
record_read2 = []
record_index1 = []
record_index2 = []
#initializing counters (and two dictionaries of counters) to count differnt kinds of reads
unknown_counter = 0
total_read_counter = 0
total_match_counter = 0
total_swapped_counter = 0
swapped_count_dict = {}
match_count_dict = {}
#looping through all input fastq files
for r1, r2, i1, i2 in zip(fh_r1, fh_r2, fh_i1, fh_i2):
    if len(record_read1) < 4:
        r1 = r1.strip()
        r2 = r1.strip()
        i1 = i1.strip()
        i2 = i2.strip()
        record_read1.append(r1)
        record_read2.append(r2)
        record_index1.append(i1)
        record_index2.append(i2)
        #only continues if all four lines of the record are appended to the list
        if len(record_read1) == 4:
            #getting the reverse comp from index 2 fastq and replacing it in the list
            record_index2[1] = revcomp(record_index2[1])
            #appending the indexes to the headers
            record_read1[0] = record_read1[0] + "-" + record_index1[1] + "-" + record_index2[1]
            record_read2[0] = record_read2[0] + "-" + record_index1[1] + "-" + record_index2[1]
            #checking if the indexes are not one the given ones, checks in Ns are in the index, and checks if the average index quality score is below the cutoff
            if (record_index1[1] not in expected_barcode) or (record_index2[1] not in expected_barcode) or ("N" in record_index1[1]) or ("N" in record_index2[1]) or (check_qscore(record_index1[3], args.cutoff)) or (check_qscore(record_index2[3], args.cutoff)):
                #writing read records to "unknown" files
                write_file(record_read1, unknown_r1)
                write_file(record_read2, unknown_r2)
                unknown_counter += 1
                total_read_counter += 1
            #if indexes are not the same (if they have been swapped)
            elif (record_index1[1] != record_index2[1]):
                #writing to swapped files
                write_file(record_read1, swap_r1)
                write_file(record_read2, swap_r2)
                #incrementing a counter in a dictionary for this pair of indexes
                up_counter(swapped_count_dict, expected_barcode[record_index1[1]] + "-" + expected_barcode[record_index2[1]])
                total_swapped_counter += 1
                total_read_counter += 1
            #if records are the same
            elif (record_index1[1] == record_index2[1]):
                #write to indexes corresponding files
                write_file(record_read1, files1[record_index1[1]])
                write_file(record_read2, files2[record_index1[1]])
                #incrementing a counter in a dictionary for this index
                up_counter(match_count_dict, expected_barcode[record_index1[1]])
                total_read_counter += 1
                total_match_counter += 1
            #empty the lists containing the records
            record_read1 = []
            record_read2 = []
            record_index1 = []
            record_index2 = []
            # For testing
            # if total_read_counter == 1000000:
            #     break
#closing all files
fh_r1.close()
fh_r2.close()
fh_i1.close()
fh_i2.close()
fh_bars.close()

for i in files1:
    files1[i].close()

for j in files2:
    files2[j].close()

#writing the report file
with open("Demux_output_report.tsv", "w") as out_fh:
    #This try and except are here since I was getting a divide by zero error. I fixed it but I'm leaving this here just in case.
    try:
        out_fh.write("Category" + "\t" + "Number of Reads" + "\t" + "Percent of Total" "\t" + "Percent of Section" + "\n")
        out_fh.write("Total Reads" + "\t" + str(total_read_counter) + "\t" + str((total_read_counter/total_read_counter)*100) + "\t" + str((total_read_counter/total_read_counter)*100) + "\n")
        out_fh.write("Matched Reads" + "\t" + str(total_match_counter) + "\t" + str((total_match_counter/total_read_counter)*100) + "\t" + str((total_match_counter/total_match_counter)*100) + "\n")
        for index in match_count_dict:
            out_fh.write(index + "\t" + str(match_count_dict[index]) + "\t" + str((match_count_dict[index]/total_read_counter)*100) + "\t" + str((match_count_dict[index]/total_match_counter)*100) + "\n")
        out_fh.write("Swapped Reads" + "\t" + str(total_swapped_counter) + "\t" + str((total_swapped_counter/total_read_counter)*100) + "\t" + str((total_swapped_counter/total_swapped_counter)*100) + "\n")
        for index in swapped_count_dict:
            out_fh.write(index + "\t" + str(swapped_count_dict[index]) + "\t" + str((swapped_count_dict[index]/total_read_counter)*100) + "\t" + str((swapped_count_dict[index]/total_swapped_counter)*100) + "\n")
        out_fh.write("Unknown Reads" + "\t" + str(unknown_counter) + "\t" + str((unknown_counter/total_read_counter)*100) + "\t" + str((unknown_counter/unknown_counter)*100) + "\n")
    except:
        print(total_match_counter)
        print(total_read_counter)
        print(total_swapped_counter)
        print(unknown_counter)
        print(swapped_count_dict)
        print(match_count_dict)

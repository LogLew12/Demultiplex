# populate this file with any fuctions you have written so far
# We suggest convert_phred() in particular :)
# Import this into your PS4.ipynb by issuing the command: import Bioinfo
def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter) - 33


def populate_list_again(file):
    """Takes fastq file and returns list of lists with individual quality scores by read by index"""
    with open(file, "r") as fh:
        ls = []
        counter = 0
        for i in range(101):
            ls.append([])
        for line in fh.readlines():
            line = line.strip()
            counter += 1
            if counter % 4 == 0:
                for index, char in enumerate(line):
                    ls[index].append(convert_phred(char))
        return ls


def gc_content(DNA):
    """Takes a string of DNA and returns an float representing the portion of the DNA that is G or C"""
    DNA = DNA.upper()
    G = DNA.count("G")
    C = DNA.count("C")
    return (G + C) / len(DNA)


def validate_base_seq(seq, RNAflag=False):
    """This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive."""
    set_DNA = set(seq.upper())
    model = {"G", "U", "C", "A", "N"} if RNAflag else {"G", "T", "C", "A", "N"}
    return set_DNA <= model


def qual_score(phred_score):
    x = 0
    for j in phred_score:
        x += convert_phred(j)
    return x / len(phred_score)


import re


def oneline_fasta(file):
    with open(file, "r") as fh:
        all_genes = {}
        for line in fh:
            line = line.strip()
            x = re.search("^>", line)
            if x:

                gene_ID = re.split("\s", line)[3]
                gene_ID = gene_ID.split(":")[1]
                protein_ID = re.split("\s", line)[0]
                protein_ID = protein_ID.split(">")[1]
                gene_name = re.search("gene_symbol:[\S]*", line)
                if gene_name == None:
                    gene_name = ""
                else:
                    gene_name = gene_name.group()
                    gene_name = gene_name.split("symbol:")[1]
                # print(gene_name)
                gene_info = (gene_ID, gene_name)
                if gene_info not in all_genes:
                    all_genes[gene_info] = [("", protein_ID)]
                else:
                    all_genes[gene_info].append(("", protein_ID))
            #     pv_was_header = True
            # elif pv_was_header == True:
            #     all_genes[gene_info].append(line)
            #     pv_was_header = False
            # elif pv_was_header == False:
            # elif len(all_genes[gene_info] == 1)
            else:
                all_genes[gene_info][len(all_genes[gene_info]) - 1] = (
                    all_genes[gene_info][len(all_genes[gene_info]) - 1][0] + line,
                    all_genes[gene_info][len(all_genes[gene_info]) - 1][1],
                )


if __name__ == "__main__":
    assert (
        validate_base_seq("AATAGAT", False) == True
    ), "Validate base seq does not work on DNA"
    assert (
        validate_base_seq("AAUAGAU", True) == True
    ), "Validate base seq does not work on RNA"
    assert validate_base_seq("TATUC", False) == False
    assert validate_base_seq("UCUGCU", False) == False
    print("Passed DNA and RNA tests")

    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")

    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
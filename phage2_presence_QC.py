import os
import copy
import sys

output = open("phage2_presence_QC.txt", "w")
script = open("../extract/single_hits.sh", "w")
rc_genomes = open("../extract/rev_complem_genomes.txt", "w")

"Species in the better category are missing less than 1463 base pairs" \
"in its alignment with the ref sequence."
better = {}

"Alignment length is exactly 14629, the same as the ref sequence," \
"and percent identity is 100%. Species here are not also in the better list."
match100percent = {}

"""Read in QC data"""
QC_file = open("../genomes-passedQC.txt", "r")
passed = QC_file.readlines()
passed[-1] = passed[-1] + "\n"

#Read in the directory that each genome is located in
locations_file = open("../alldata-asmstats_checkm.tsv", "r")
locations_file.readline()
line0 = locations_file.readline()
locations = {}
while line0 != "":
    genome = line0.split()[0]
    dir0 = line0.split()[-1]
    locations[genome] = dir0
    line0 = locations_file.readline()

def do(setting):
    full_path = "/gsap/archive-bacterial/Projects/EnteroGenome/phage2/BLAST_results/" + setting
    for f in os.listdir(full_path):
        sample = ""
        if ".txt" in f:
            sample = f[0:f.index(".txt")]
        if not os.path.isdir(full_path + f) and sample + "\n" in passed:
            blast_file = open(full_path + f, "r")
            for i in range(0, 5):
                blast_file.readline()

            maxlen = 0
            maxiden = 0
            qstarts = []
            qends = []
            sstarts = []
            s_ends = []
            contigs = []
            all_short = True
            calculate_coverage = True
            line0 = blast_file.readline()
            while line0 != "":
                if not line0.startswith("#"):
                    length = int(line0.split()[3])
                    iden = float(line0.split()[2])
                    if length >= 2000:
                        all_short = False
                    if length >= 14629 and iden > 99.99:
                        match100percent[sample] = (line0.split()[1], line0.split()[8], line0.split()[9])
                        calculate_coverage = False
                        if sample in better:
                            print(sample, "goes into the better and 100% categories depending on the setting")
                            print(setting)
                        break
                    elif length >= 14629:
                        if sample in match100percent:
                            print(sample, "goes into the better and 100% categories depending on the setting")
                            print(setting)
                        if sample in better:
                            if len(better[sample]) == 5 and iden > better[sample][1]:
                                better[sample] = (length, iden, line0.split()[1], line0.split()[8], line0.split()[9])
                            elif len(better[sample]) == 2:
                                print(sample, "has One or Multiple hits in Better depending on the setting")
                        else:
                            better[sample] = (length, iden, line0.split()[1], line0.split()[8], line0.split()[9])
                        calculate_coverage = False
                        break

                    qstarts.append(int(line0.split()[6]))
                    qends.append(int(line0.split()[7]))
                    sstarts.append(line0.split()[8])
                    s_ends.append(line0.split()[9])
                    contigs.append(line0.split()[1])
                    if length > maxlen:
                        maxlen = length
                        maxiden = iden
                line0 = blast_file.readline()
            if not all_short and calculate_coverage and maxiden >= 80:
                "The r stands for remainder fragments. Each element is a" \
                "[start, end] pair of Phage2 sections that haven't been covered yet."
                r = [[1, 14629]]
                coordinates = []
                for i in range(len(qstarts)):
                    if qends[i] - qstarts[i] > 400:
                        new_r = copy.deepcopy(r)
                        for fragment in r:
                            """If the hit is contained within a remainder fragment"""
                            if qstarts[i] >= fragment[0] and qends[i] <= fragment[1]:
                                if qstarts[i] != fragment[0]:
                                    new_r.append([fragment[0], qstarts[i]-1])
                                if qends[i] != fragment[1]:
                                    new_r.append([qends[i] + 1, fragment[1]])
                                new_r.remove(fragment)
                                coordinates.append(contigs[i] + ": " + str(sstarts[i]) + ", " + str(s_ends[i]) + "\t")
                            elif qstarts[i] >= fragment[0] and qstarts[i] <= fragment[1]:
                                """If the hit starts within a remainder fragment and ends outside it"""
                                """The following line checks if the new fragment you're adding is invalid ([2, 1] for example)"""
                                if qstarts[i] != fragment[0]:
                                    new_r.append([fragment[0], qstarts[i]-1])
                                new_r.remove(fragment)
                                coordinates.append(contigs[i] + ": " + str(sstarts[i]) + ", " + str(s_ends[i]) + "\t")
                            elif qends[i] >= fragment[0] and qends[i] <= fragment[1]:
                                """If the hit starts outside remainder fragments and ends inside one"""
                                """The following line checks if the new fragment you're adding is invalid ([2, 1] for example)"""
                                if qends[i] != fragment[1]:
                                    new_r.append([qends[i] + 1, fragment[1]])
                                new_r.remove(fragment)
                                coordinates.append(contigs[i] + ": " + str(sstarts[i]) + ", " + str(s_ends[i]) + "\t")
                            elif qstarts[i] <= fragment[0] and qends[i] >= fragment[1]:
                                """If the remainder fragment is contained within the hit"""
                                new_r.remove(fragment)
                                coordinates.append(contigs[i] + ": " + str(sstarts[i]) + ", " + str(s_ends[i]) + "\t")
                    r = new_r
                    if r == []:
                        break
                """Now determine the coverage by adding up the remainder fragments"""
                inverse_coverage = 0
                for fragment in r:
                    inverse_coverage += fragment[1] - fragment[0] + 1
                if inverse_coverage < 1463:
                    if sample in match100percent:
                        print(sample, "goes into the better and 100% categories depending on the setting")
                        print(setting)
                    if sample in better:
                        if len(better[sample]) == 2 and inverse_coverage < better[sample][0]:
                            better[sample] = (inverse_coverage, coordinates)
                        elif len(better[sample]) == 5:
                            print(sample, "has One or Multiple hits in Better depending on the setting")
                    else:
                        better[sample] = (inverse_coverage, coordinates)



do("blastn/full_match/")
do("megablast/full_match/")
do("dc-megablast/full_match/")
do("blastn/")
do("megablast/")
do("dc-megablast/")

output.write("100% matches:\n(Alignment length is exactly 14629, the same as the Phage 2 reference sequence, "
             "and percent identity is 100%. There is no overlap between this list and the next list.) The columns are: genome, contig, start coordinates, and end coordinates.\n")
for sample in match100percent:
    output.write(sample + "\t" + match100percent[sample][0] + "\t" + str(match100percent[sample][1]) + "\t" + str(match100percent[sample][2]) + "\n")
    """Writing to the extract script"""
    if int(match100percent[sample][1]) <= int(match100percent[sample][2]):
        phage2_start = match100percent[sample][1]
        phage2_end = match100percent[sample][2]
        left_start = int(phage2_start) - 11378
        if left_start < 1:
            left_start = 1
        left_start = str(left_start)
        right_end = str(int(phage2_end) + 7893)
    else:
        phage2_start = match100percent[sample][2]
        phage2_end = match100percent[sample][1]
        left_start = int(phage2_start) - 7893
        if left_start < 1:
            left_start = 1
        left_start = str(left_start)
        right_end = str(int(phage2_end) + 11378)
        rc_genomes.write(sample + "\n")
    dir0 = locations[sample]
    if dir0 == "ncbi":
        script.write(
            "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/ncbi/all-genomes/genome/"
            + sample + "_genomic.fna " + match100percent[sample][0] + ":" + left_start + "-" + phage2_start + " > " +
            sample + "_left.fa \n")
        script.write(
            "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/ncbi/all-genomes/genome/"
            + sample + "_genomic.fna " + match100percent[sample][0] + ":" + phage2_start + "-" + phage2_end + " > " +
            sample + ".fa \n")
        script.write(
            "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/ncbi/all-genomes/genome/"
            + sample + "_genomic.fna " + match100percent[sample][0] + ":" + phage2_end + "-" + right_end + " > " +
            sample + "_right.fa \n")
    elif dir0 == "divlibrary" or dir0 == "zoo2/divlibrary" or dir0 == "zoo2":
        script.write(
            "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/divlibrary/genome/"
            + sample + ".genome.fa " + match100percent[sample][0] + ":" + left_start + "-" + phage2_start + " > " +
            sample + "_left.fa \n")
        script.write(
            "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/divlibrary/genome/"
            + sample + ".genome.fa " + match100percent[sample][0] + ":" + phage2_start + "-" + phage2_end + " > " +
            sample + ".fa \n")
        script.write(
            "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/divlibrary/genome/"
            + sample + ".genome.fa " + match100percent[sample][0] + ":" + phage2_end + "-" + right_end + " > " +
            sample + "_right.fa \n")
    elif dir0 == "zoo3":
        script.write(
            "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/zoo3/annotation/prokka/"
            + sample + "/" + sample + ".fna " + match100percent[sample][0] + ":" + left_start + "-" + phage2_start + " > " +
            sample + "_left.fa \n")
        script.write(
            "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/zoo3/annotation/prokka/"
            + sample + "/" + sample + ".fna " + match100percent[sample][0] + ":" + phage2_start + "-" + phage2_end + " > " +
            sample + ".fa \n")
        script.write(
            "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/zoo3/annotation/prokka/"
            + sample + "/" + sample + ".fna " + match100percent[sample][0] + ":" + phage2_end + "-" + right_end + " > " +
            sample + "_right.fa \n")
    else:
        print(sample + " not found in file locations excel file")
        sys.exit()

output.write("Each hit was subtracted off of the reference sequence, and the sum of the remaining reference fragments "
             "calculated. If this value is less than 1463, which is 10% of the length of the reference sequence, the "
             "genome is present in this list. In addition, the longest hit must have percent identity greater than 80%. The columns are: contig, start coordinates, and end coordinates for Single Hit genomes, and a list of contig coordinate pairs for Multiple Hit genomes.\n")
for sample in better:
    if len(better[sample]) == 5:
        output.write(sample + "\tSingle hit\t" + str(better[sample][2]) + "\t" + str(better[sample][3]) + "\t" +
                     str(better[sample][4]) + "\n")
        """Writing to the extract script"""
        if int(better[sample][3]) <= int(better[sample][4]):
            phage2_start = better[sample][3]
            phage2_end = better[sample][4]
            left_start = int(phage2_start) - 11378
            if left_start < 1:
                left_start = 1
            left_start = str(left_start)
            right_end = str(int(phage2_end) + 7893)
        else:
            phage2_start = better[sample][4]
            phage2_end = better[sample][3]
            left_start = int(phage2_start) - 7893
            if left_start < 1:
                left_start = 1
            left_start = str(left_start)
            right_end = str(int(phage2_end) + 11378)
            rc_genomes.write(sample + "\n")
        dir0 = locations[sample]
        if dir0 == "ncbi":
            script.write(
                "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/ncbi/all-genomes/genome/"
                + sample + "_genomic.fna " + better[sample][2] + ":" + left_start + "-" + phage2_start + " > " +
                sample + "_left.fa \n")
            script.write(
                "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/ncbi/all-genomes/genome/"
                + sample + "_genomic.fna " + better[sample][2] + ":" + phage2_start + "-" + phage2_end + " > " +
                sample + ".fa \n")
            script.write(
                "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/ncbi/all-genomes/genome/"
                + sample + "_genomic.fna " + better[sample][2] + ":" + phage2_end + "-" + right_end + " > " +
                sample + "_right.fa \n")
        elif dir0 == "divlibrary" or dir0 == "zoo2/divlibrary" or dir0 == "zoo2":
            script.write(
                "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/divlibrary/genome/"
                + sample + ".genome.fa " + better[sample][2] + ":" + left_start + "-" + phage2_start + " > " +
                sample + "_left.fa \n")
            script.write(
                "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/divlibrary/genome/"
                + sample + ".genome.fa " + better[sample][2] + ":" + phage2_start + "-" + phage2_end + " > " +
                sample + ".fa \n")
            script.write(
                "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/divlibrary/genome/"
                + sample + ".genome.fa " + better[sample][2] + ":" + phage2_end + "-" + right_end + " > " +
                sample + "_right.fa \n")
        elif dir0 == "zoo3":
            script.write(
                "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/zoo3/annotation/prokka/"
                + sample + "/" + sample + ".fna " + better[sample][2] + ":" + left_start + "-" + phage2_start + " > " +
                sample + "_left.fa \n")
            script.write(
                "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/zoo3/annotation/prokka/"
                + sample + "/" + sample + ".fna " + better[sample][2] + ":" + phage2_start + "-" + phage2_end + " > " +
                sample + ".fa \n")
            script.write(
                "samtools faidx /gsap/archive-bacterial/Projects/EnteroGenome/zoo3/May2021/data/zoo3/annotation/prokka/"
                + sample + "/" + sample + ".fna " + better[sample][2] + ":" + phage2_end + "-" + right_end + " > " +
                sample + "_right.fa \n")
        else:
            print(sample + " not found in file locations excel file")
            sys.exit()
    else:
        output.write(sample + "\tMultiple hits\t")
        for i in better[sample][1]:
            output.write(i)
        output.write("\n")

import sys

"""file0 = open("alignment2.afa", "r")
output = open("alignment_subset.txt", "w")
line = file0.readline()
while line != "":
	line = file0.readline()
	output.write(line[0:100])
	while not line.startswith(">") and line != "":
		line = file0.readline()"""

"""for i in range(10):
	print(i)
	file1 = open("alignment_subset.txt", "r")
	first_line = file1.readline()
	line = first_line
	base = first_line[i]
	while line != "":
		if line[i] != base:
			print(first_line)
			print(line)
			print()
		line = file1.readline()"""

#To get the start coordinate, find the coordinate in the "even_better" annotations file, subtract one, and add the number of dashes by the beginning of the nth gene calculated below. For the end coordinate, don't subtract one, then add the number of dashes by the end of the nth gene
#Print first gene
file0 = open("alignment2.afa", "r")
line = file0.readline()
while "006494835" not in line:
	line = file0.readline()
line = file0.readline()
refseq = ""
while not line.startswith(">"):
	refseq += line.replace("\n", "")
	line = file0.readline()
first_gene = refseq[0:476]
#print(first_gene)
gene2 = refseq[539:1024] #The 541: -1 for 1-indexing, +19 for the number of dashes before the stop codon of the 2nd gene. The 1024: +21 for the number of dashes before the start codon of the 2nd gene
#print(gene2)
gene3 = refseq[1248:1384]
#print(gene3)
#print("Number of dashes by the beginning of the 2nd gene including intergenic regions:", refseq[0:522].count("-"))
#print("Number of dashes by the end of the 2nd gene including intergenic regions:", refseq[0:1005].count("-"))
#print("Number of dashes by the beginning of the 3rd gene including intergenic regions:", refseq[0:1248].count("-"))
#print("Number of dashes by the end of the 3rd gene including intergenic regions:", refseq[0:1384].count("-"))
gene = refseq[1414:2199]
print(gene)
#print("Number of dashes by the beginning of the nth gene including intergenic regions:", refseq[0:13699].count("-"))
#print("Number of dashes by the end of the nth gene including intergenic regions:", refseq[0:14798].count("-"))
sys.exit()

"""codon_lookup = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L", "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "UAU":"Y", "UAC":"Y", "UAA":"stop", "UAG":"stop", "UGU":"C", "UGC":"C", "UGA":"stop", "UGG":"W", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P", "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M", "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T", "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R", "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A", "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

def complement(base):
	if base == "A":
		return "U"
	elif base == "T":
		return "A"
	elif base == "C":
		return "G"
	elif base == "G":
		return "C"
	elif base == "-":
		return "-"

genomes = []
proteins = []
file0 = open("alignment2.afa", "r")
line = file0.readline()
#while "006494835" not in line:
#	line = file0.readline()
while line != "":
	genomes.append(line[1:-1])
	sequence = ""
	protein = ""
	matrix_row = ""
	line = file0.readline()
	while not line.startswith(">") and line != "":
		sequence += line.replace("\n", "")
		line = file0.readline()
	first_gene = sequence[0:476]  #This is where AUG starts
	#if sequence[230] != "-":
	#	print(sequence[475:229:-1])
	for i in range(476-3, -1, -3):
		codon = first_gene[i : i+3]
		#print(codon)
		#sys.exit()
		anticodon = complement(codon[2]) + complement(codon[1]) + complement(codon[0])
		if "-" in anticodon:
			protein += "-"
		else:
			protein += codon_lookup[anticodon]
	#print(protein)
	proteins.append(protein)
"""
#It's 2774 for both!
#print(len(genomes))
#print(len(proteins))

#ref_index = genomes.index("GCA_006494835.1_ASM649483v1")
#ref_protein = proteins[ref_index]
#nonsyn_matrix = []

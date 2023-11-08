import sys

#Skips any codons with N's
#All the nonsynonymous mutations (but not including deletions, and insertions aren't handled in this script), their base pair position will be reported as the position of the start of the codon. So if the ref sequence is AAACCC, then for AAATCC, AAACTC, and AAACCT this script will say there's a mutation at position 3 (0-indexed)

codon_lookup = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L", "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "UAU":"Y", "UAC":"Y", "UAA":"stop", "UAG":"stop", "UGU":"C", "UGC":"C", "UGA":"stop", "UGG":"W", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P", "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M", "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T", "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R", "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A", "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

def complement(base):
	if base == "A":
		return "U"
	elif base == "T":
		return "A"
	elif base == "C":
		return "G"
	elif base == "G":
		return "C"

for gene_num in range(1276, 1294):
	#Translate and compare with the reference's protein sequence, this will find nonsynonymous mutations. Within this you can count how many are stop codons
	#Have to read in insertions files, stop translating after the first insertion in each marked genome, remember the exceptions for the insertions in EF1288 that don't cause a frameshift
	#Deletions, or dashes in non-ref genomes. Will error with codon_lookup, just save it and say it's a frameshift, unless there's a multiple of 3 dashes, then print an alert for that
	gene_alignment = open("alignment_no_insertions_EF" + str(gene_num) + ".afa", "r")
	insertions_file = open("insertions_EF" + str(gene_num) + ".txt", "r")
	output = open("/gsap/archive-bacterial/Projects/EnteroGenome/phage2/all_faecalis_whole_genome_phylogeny/nonsynonymous/EF" + str(gene_num) + ".txt", "w")
	line = insertions_file.readline()
	insertions = {}
	if "No insertions found!" not in line:
		while line != "":
			position = int(line.split()[1][0:-1])
			genome_list = insertions_file.readline().split()
			if not (gene_num == 1288 and genome_list[0] in ["GCA_905123795.1_28157_4_332", "GCA_902161785.1_25426_7_322", "GCA_905122585.1_28099_2_18"]):
				insertions[position] = genome_list
			line = insertions_file.readline()
			line = insertions_file.readline()

	#Save a 2D matrix of the alignment. Also have a simple list of the genomes, and they'll be in the same order as the rows of the matrix.
	genomes = []
	dna_seqs = []
	line = gene_alignment.readline()
	while line != "":
		if not line.startswith(">"):
			print("Error")
			sys.exit()
		genomes.append(line[1:-1])
		sequence = gene_alignment.readline().replace("\n", "")
		dna_seqs.append(sequence)
		line = gene_alignment.readline()
	num_genomes = len(genomes)
	ref_index = genomes.index("GCA_006494835.1_ASM649483v1")
	ref_dna = dna_seqs[ref_index]

	#Translate the reference sequence
	#The first 2 genes are in the reverse orientation. In the alignment_no_insertions files, I reversed them so I can still read them left to right here, but I still have to complement the sequence
	ref_protein = ""	
	check = False
	for i in range(0, len(ref_dna), 3):
		codon = ref_dna[i : i+3]
		if gene_num == 1276 or gene_num == 1277:
			anticodon = complement(codon[0]) + complement(codon[1]) + complement(codon[2])
		else:
			if codon[0] == "T":
				anticodon = "U"
			else:
				anticodon = codon[0]
			if codon[1] == "T":
				anticodon += "U"
			else:
				anticodon += codon[1]
			if codon[2] == "T":
				anticodon += "U"
			else:
				anticodon += codon[2]
		ref_protein += codon_lookup[anticodon]
		if i == len(ref_dna) - 3:
			check = True
	if not check:
		print("Error: not all of the ref sequence has been translated")
		sys.exit()
	
	#if gene_num == 1291:
	#	for pos in [1545, 1640, 1542, 1556, 1561, 1574, 1579, 1581, 1589, 1592, 1595, 1600, 1601, 1602, 1604, 1607, 1614, 1617, 1625, 1629, 1631, 1643, 1652, 1663, 1664, 1692, 1659, 1701, 1544, 1586, 1587, 1597, 1616, 1541, 1705, 1583, 1598, 1548, 1599, 1637, 1676, 1673, 1709, 1699, 1549, 1563, 1565, 1571, 1577, 1620, 1622, 1623, 1627, 1633, 1634, 1638, 1647, 1660, 1665, 1688, 1691, 1703, 1707, 1630, 1646, 1649, 1700]:
	#		print(pos//3)
	#		print("Position", pos, "is", ref_protein[pos//3])

	#Translate the rest of the genomes
	done = []
	for i in range(num_genomes):
		genome = genomes[i]
		seq = dna_seqs[i]
		stopping_point = len(seq)
		#print(len(ref_dna))
		#print(stopping_point)
		for j in insertions:
			if genome in insertions[j]:
				stopping_point = j
				if stopping_point % 3 != 0:
					stopping_point -= 3
				break
		ref_protein_pos = 0
		mutations = []
		for j in range(0, stopping_point, 3):
			codon = seq[j : j+3]
			if "N" in codon:
				dummy = 0
			elif "-" in codon:
				del_pos = j + codon.index("-")
				if seq[del_pos] != "-":
					print("Error on line 94, EF" + str(gene_num), "genome", genome)
				k = del_pos
				while seq[k] == "-":
					k += 1
				if (k-del_pos) % 3 == 0:
					output.write(genome + "\tNo frameshift deletion. " + str(k-del_pos) + " bp deletion at position " + str(del_pos) + "\n")
				else:
					output.write(genome + "\tFrameshift deletion. " + str(k-del_pos) + " bp deletion at position " + str(del_pos) + "\n")
					break
				#else:
				#	print("Non-frameshift-causing deletion found. Gene: EF" + str(gene_num), "Genome:", genome, "Position:", del_pos, "-", k)
			else:
				if gene_num == 1276 or gene_num == 1277:
					anticodon = complement(codon[0]) + complement(codon[1]) + complement(codon[2])
					#print(anticodon)
				else:
					if codon[0] == "T":
						anticodon = "U"
					else:
						anticodon = codon[0]
					if codon[1] == "T":
						anticodon += "U"
					else:
						anticodon += codon[1]
					if codon[2] == "T":
						anticodon += "U"
					else:
						anticodon += codon[2]
				amino_acid = codon_lookup[anticodon]
				#print(amino_acid)
				ref_amino_acid = ref_protein[ref_protein_pos]
				#print(ref_amino_acid)
				if ref_amino_acid == "s" and amino_acid != "stop":
					output.write(genome + "\t" + "nonstop mutation\n")
					#print("ref is stop but genome isn't")
					#print("EF" + str(gene_num), genome)
					#print("Reference seq:", ref_dna)
					#print("The genome:", seq)
					#sys.exit()
					break
				if amino_acid == "stop" and ref_amino_acid != "s":
					output.write(genome + "\t" + "nonsense mutation at position " + str(j) + "-" + str(j+2) + "\n")
					#print("Nonsense")
					#print("EF" + str(gene_num), genome, str(j) + "-" + str(j+2))
					#print("Reference seq:", ref_dna)
					#print("The genome:", seq)
					#sys.exit()
					break
				if ref_amino_acid == "s" and amino_acid == "stop":
					#print("It works")
					break
				if ref_amino_acid != amino_acid:
					mutations.append(str(j))
					if j not in done:
						if gene_num == 1291:
							print("Nonsyn mutation at position", j, "ref amino acid is", ref_amino_acid)
							done.append(j)
					"""if len(mutations) == 2:
						print(mutations)
						print("EF" + str(gene_num), genome)
						print("Reference seq:", ref_dna)
						print("The genome:", seq)
						sys.exit()"""
			ref_protein_pos += 1
		
		if len(mutations) > 0:
			output.write(genome + "\t" + "locations of nonsynonymous mutations, not including nonsense: " + ",".join(mutations) + "\n")

		"""min_len = min(len(ref_protein), len(protein))
		j = 0
		while j < min_len:
			residue1 = ref_protein[j]
			residue2 = protein[j]
			if residue1 == "s" and residue2 != "s":
				output.write(genome + "\t" + "ref stop codon mutated to amino acid\n")
				break
			if residue2 == "s" and residue1 != "s":
				dna_position = 3*j
				output.write(genome + "\t" + "nonsense mutation at position " + str(dna_position) + "-" + str(dna_position+2) + "\n")
				break
			if residue1 != residue2:
				output.write(genome + "\t" + residue1 + " at position " + str(j) + " mutated to " + residue2 + "\n")
			if residue1 == "s":
				break
			j += 1"""

"""proteins = []
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

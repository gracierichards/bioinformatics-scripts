import sys

#Save a 2D matrix of the alignment. Also have a simple list of the genomes, and they'll be in the same order as the rows of the matrix.
genomes = []
dna_seqs = []

alignment = open("../all_faecalis_whole_genome_phylogeny/snps/rep_alignment.afa", "r")
line = alignment.readline()
while line != "":
	if not line.startswith(">"):
		print("Error")
		sys.exit()
	genomes.append(line[1:-1])
	sequence = alignment.readline().replace("\n", "")
	dna_seqs.append(sequence)
	line = alignment.readline()

num_genomes = len(genomes)
ref_index = genomes.index("GCA_006494835.1_ASM649483v1")
ref_sequence = dna_seqs[ref_index]

#end_index is the length of the gene, not len-1.
def find_insertions(gene_name, start_index, end_index, reverse_orient):
	seqs_without_insertions = [""] * num_genomes
	
	alignment_without_insertions = open("alignment_no_insertions_" + gene_name + ".afa", "w")
	insertions_file = open("insertions_" + gene_name + ".txt", "w")

	#The keys are the indices of the insertions, the values are a list of genomes that have an insertion there
	#The index is relative to the beginning of the gene whether the gene is in forward or reverse orientation. It's 0-indexed
	insertions = {}
	
	if reverse_orient:
		position = 0
		if start_index == 0:
			gene = ref_sequence[end_index-1:0:-1] + ref_sequence[0]
		else:
			gene = ref_sequence[end_index-1:start_index-1:-1]
	else:
		position = start_index
		gene = ref_sequence[start_index:end_index]
		#TEMP CODE FOR FINDING EF1291 POSITIONS
		positions2 = [12903, 13005, 12900, 12914, 12919, 12932, 12937, 12939, 12954, 12957, 12960, 12965, 12966, 12967, 12969, 12972, 12979, 12982, 12990, 12994, 12996, 13008, 13017, 13028, 13029, 13057, 13024, 13066, 12902, 12951, 12952, 12962, 12981, 12899, 13070, 12948, 12963, 12906, 12964, 13002, 13041, 13038, 13074, 13064, 12907, 12921, 12923, 12929, 12935, 12985, 12987, 12988, 12992, 12998, 12999, 13003, 13012, 13025, 13030, 13053, 13056, 13068, 13072, 12995, 13011, 13014, 13065]
		positions3 = [12903, 13005, 12900, 12914, 12919, 12932, 12937, 12939, 12954, 12957, 12960, 12965, 12966, 12967, 12969, 12972, 12979, 12982, 12990, 12994, 12996, 13008, 13017, 13028, 13029, 13057, 13024, 13066, 12902, 12951, 12952, 12962, 12981, 12899, 13070, 12948, 12963, 12906, 12964, 13002, 13041, 13038, 13074, 13064, 12907, 12921, 12923, 12929, 12935, 12985, 12987, 12988, 12992, 12998, 12999, 13003, 13012, 13025, 13030, 13053, 13056, 13068, 13072, 12995, 13011, 13014, 13065]
		
		#for i in positions2:
		#	print("Position", str(i) + ":", ref_sequence[i])
		
	#Looking for insertions and deleting gapped columns during this loop
	#As it writes the first gene to seqs_without_insertions, the sequences are reversed
	while "-" in gene:
		dash_index = gene.index("-")
		
		#Testing the coordinates of the subsequence being written to file for forward orientation
		#if not reverse_orient:
		#	print("Index of dash is", position+dash_index)

		position_after_dash = end_index-1-position
		if not reverse_orient:
			insertions[position-start_index-1] = []
			#Insertions relative to the start of the whole Phage 2 sequence:
			#insertions[position+dash_index] = []
			#Testing if all dashes have been written to insertions
			#print(position+dash_index)

		#Testing if position_after_dash+1 is actually a dash in the reference sequence
		#if reverse_orient and position != 0:
			#print(ref_sequence[position_after_dash+1])
			#print(position_after_dash+1)
		#if not reverse_orient:
		#	print(ref_sequence[position + dash_index])
		#once = True
		for i in range(num_genomes):
			if reverse_orient:
				#Testing if all the dashes have been written to insertions
				#if position != 0 and once:
				#	print(position_after_dash+1, "has been checked")
				#	once = False

				if position != 0 and dna_seqs[i][position_after_dash+1] != "-":
					if dna_seqs[i][position_after_dash+1] in ["A", "C", "G", "T"]:
						insertions[position - 1].append(genomes[i])
						#Insertions relative to the start of the whole Phage 2 sequence:
						#insertions[position_after_dash+1].append(genomes[i])
					else:
						print("Position", position_after_dash+1, "is a dash in the reference sequence, but in genome", genomes[i], "it's neither -, A, C, G, or T")
				seqs_without_insertions[i] += dna_seqs[i][position_after_dash: position_after_dash-dash_index : -1]
			else:
				#Testing if all the dashes have been written to insertions
				#if once:
				#	print(position+dash_index, "has been checked")
				#	once = False

				if dna_seqs[i][position+dash_index] != "-":
					if dna_seqs[i][position+dash_index] in ["A", "C", "G", "T"]:
						insertions[position - start_index - 1].append(genomes[i])
						#Insertions relative to the start of the whole Phage 2 sequence:
						#insertions[position+dash_index].append(genomes[i])
					else:
						print("Position", position+dash_index, "is a dash in the reference sequence, but in genome", genomes[i], "it's neither -, A, C, G, or T")
				seqs_without_insertions[i] += dna_seqs[i][position:position+dash_index]

		gene = gene[dash_index + 1:]
		position += dash_index + 1
		
		#TEMP CODE FOR FINDING EF1291 POSITIONS		
		for i in range(len(positions2)):
			if positions2[i] > position:
				positions3[i] -= 1

		if reverse_orient:
			insertions[position - 1] = []
			#Insertions relative to the start of the whole Phage 2 sequence:
			#insertions[end_index-position] = []
	
	#Doing the same steps for the sequence after the last dash, as the loop terminates when it sees the last dash
	position_after_dash = end_index-1-position
	#Testing if all the dashes have been written to insertions
	#if reverse_orient:
	#	print(position_after_dash+1)
	#once = True
	for i in range(num_genomes):
		if reverse_orient:
			#Testing if all the dashes have been written to insertions
			#if once:
			#	print(position_after_dash+1, "has been checked")
			#	once = False

			if dna_seqs[i][position_after_dash+1] != "-":
				if dna_seqs[i][position_after_dash+1] in ["A", "C", "G", "T"]:
					insertions[position - 1].append(genomes[i])
					#Insertions relative to the start of the whole Phage 2 sequence:
					#insertions[position_after_dash+1].append(genomes[i])
				else:
					print("Position", position_after_dash+1, "is a dash in the reference sequence, but in genome", genomes[i], "it's neither -, A, C, G, or T")
			if start_index == 0:
				seqs_without_insertions[i] += dna_seqs[i][position_after_dash:0:-1] + dna_seqs[i][0]
			else:
				seqs_without_insertions[i] += dna_seqs[i][position_after_dash: start_index-1 : -1]
		else:
			seqs_without_insertions[i] += dna_seqs[i][position:end_index]
	
	#TEMP CODE TO FIND EF1291 POSITIONS
	#print(seqs_without_insertions[ref_index][12944-start_index-21:12949-start_index-21])
	#sys.exit()
	#print("12900 after removing dashes is", positions2[2])
	#print(seqs_without_insertions[ref_index][positions2[2]-start_index-4:positions2[2]-start_index+5])
	print("BREAK")
	#print(seqs_without_insertions[ref_index][positions2[0] - start_index - 4: positions2[2] - start_index + 5])
	for i in positions3:
		print("Position", str(i) + ":", seqs_without_insertions[ref_index][i-start_index])
	print([i - start_index for i in positions3])
	"""print("Final position is", positions2[1])
	print("Length of ref seq is", end_index - start_index)
	print("Length of ref seq without insertions is", len(seqs_without_insertions[ref_index]))
	print("Indexing position", positions2[1] - start_index, "of ref seq without insertions")
	print(seqs_without_insertions[ref_index][positions2[1] - start_index])
	print(ref_sequence[start_index:end_index])
	print(ref_sequence[start_index:end_index].replace("-", ""))
	print(seqs_without_insertions[ref_index])"""


	#Writing to files
	for i in range(num_genomes):
		alignment_without_insertions.write(">" + genomes[i] + "\n")
		alignment_without_insertions.write(seqs_without_insertions[i] + "\n")
	all_zero = True
	for i in insertions:
		if len(insertions[i]) != 0:
			insertions_file.write("Position " + str(i) + ":\n")
			for genome in insertions[i]:
				insertions_file.write(genome + " ")
			insertions_file.write("\n\n")
			all_zero = False
	if all_zero:
		insertions_file.write("No insertions found!")

find_insertions("EF1276", 0, 476, True)
find_insertions("EF1277", 539, 1024, True)
find_insertions("EF1278", 1248, 1384, False)
find_insertions("EF1279", 1414, 2199, False)
find_insertions("EF1280", 2217, 3066, False)
find_insertions("EF1281", 3068, 3161, False)
find_insertions("EF1282", 3153, 3545, False)
find_insertions("EF1283", 3564, 3977, False)
find_insertions("EF1284", 4225, 4629, False)
find_insertions("EF1285", 4641, 5154, False)
find_insertions("EF1286", 5188, 5550, False)
find_insertions("EF1287", 5567, 5937, False)
find_insertions("EF1288", 5925, 8905, False)
find_insertions("EF1289", 8906, 9833, False)
find_insertions("EF1290", 9849, 11312, False)
find_insertions("EF1291", 11344, 13135, False)
find_insertions("EF1292", 13158, 13554, False)
find_insertions("EF1293", 13699, 14798, False)

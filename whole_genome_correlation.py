import os
import sys
from scipy import stats
#from pandas import DataFrame
#import seaborn as sns
#from seaborn import clustermap
#from seaborn import histplot
#import matplotlib.pyplot as plt

gene_coords = [[0, 475], [539, 1023], [1248, 1383], [1414, 2198], [2217, 3065], [3068, 3160], [3153, 3544], [3564, 3976], [4225, 4628], [4641, 5153], [5188, 5549], [5567, 5936], [5925, 8904], [8906, 9832], [9849, 11311], [11344, 13134], [13158, 13553], [13699, 14797]]
epa_orthogroups = {"021156":"EF2163", "025142":"EF2164", "016357":"EF2165", "014989":"EF2166", "021745":"EF2167", "022911":"EF2168", "022421":"EF2169", "011034":"EF2170", "021466":"EF2171", "020215":"EF2172", "019203":"EF2173", "016915":"EF2174", "025468":"EF2175", "025486":"EF2176", "014339":"EF2177", "015147":"EF2178", "021778":"EF2179", "024179":"EF2180", "015153":"EF2181", "014357":"EF2182", "016188":"EF2183", "024272":"EF2184", "019203":"EF2185", "012218":"EF2186", "016704":"EF2187", "019246":"EF2188", "021095":"EF2189", "014264":"EF2190", "023040":"EF2191", "014569":"EF2192", "022484":"EF2193", "017780":"EF2194", "015140":"EF2195", "018237":"EF2196", "014550":"EF2197", "017253":"EF2198", "017281":"EF2199"}
gene_functions = {}

gene_cols = {}
synerclust_file = open("results400/results/cluster_dist_per_genome.txt", "r")
#output2 = open("orthogroups_found_in_all_genomes.txt", "w")
#output3 = open("genome_wide_correlations.txt", "w")
#output4 = open("all_genes_with_EF1291_correlations.txt", "w")
#output3.write("Orthogroup\tAnnotation for orthogroup\tPhage 2 gene\tCorrelation\tepa gene\n")
#output4.write("Orthogroup\tAnnotation for orthogroup\tPhage 2 gene\tCorrelation\tepa gene\n")
header = synerclust_file.readline()
genomes_in_order = header.split()[2:]
i = genomes_in_order.index("GCA_000007785.1_ASM778v1.fna.ref")
line = synerclust_file.readline()
while line != "":
	digits = line.split()
	cluster = digits[0]
	column = []
	for j in range(2, len(digits)):
		#Deleting the ref genome column because it's not in the SNP file
		if j-2 != i:
			digit = digits[j]
			if digit == "0":
				column.append(0)
			elif digit == "1":
				column.append(1)
			else:
				column.append(1)
	#If the whole column has the same value, omit this orthogroup
	all_same = True
	val = column[0]
	for j in column:
		if j != val:
			all_same = False
	if all_same:
		x = 0
		#if val == 0:
		#	output2.write(cluster + " is absent from all genomes except for GCA_000007785.1_ASM778v1\n")
		#else:
		#	output2.write(cluster + "\n")
	else:
		gene_cols[cluster] = column
	line = synerclust_file.readline()

snp_annot_cols = {}
all_data = {}
snp_file = open("../../../all_faecalis_whole_genome_phylogeny/snps/sorted_150bp_binary_heatmap.txt", "r")
line = snp_file.readline()
while line != "DATA\n":
	line = snp_file.readline()
	if line.startswith("FIELD_LABELS"):
		genes_in_order = line.split()[1:]
for gene in genes_in_order:
	snp_annot_cols[gene] = []
line = snp_file.readline()
while line != "":
	genome = line.split()[0].replace("'", "")
	genome2 = genome[0:-4]
	all_data[genome2] = line.split()[1:]
	line = snp_file.readline()
for item in genomes_in_order:
	if item != "GCA_000007785.1_ASM778v1.fna.ref":
		#Remove all .fna and .fna.ref
		genome = item[0 : item.index(".fna")]
		#if genome not in all_data:
		#	print(genomes_in_order.index(item))
		#	print("Number of genomes =", len(genomes_in_order))
		for i in range(len(genes_in_order)):
			gene = genes_in_order[i]
			digit = int(all_data[genome][i])
			snp_annot_cols[gene].append(digit)
row_names = genes_in_order

#pearsonr will error if every value in a column is the same, so this part checks for that and prints the columns that are constant. They will also be omitted from the heatmap
phage2_genes_to_skip = []
for phage2_gene in snp_annot_cols:
	column = snp_annot_cols[phage2_gene]
	all_same = True
	val = column[0]
	for i in column:
		if i != val:
			all_same = False
			break
	if all_same:
		print("The matrix column,", phage2_gene, "consists of all", str(val) + "'s. Omitting from heatmap.")
		phage2_genes_to_skip.append(phage2_gene)

all_correlations = []
EF1291_correlations = []
for cluster in gene_cols:
	column = []
	printed = False
	for chunk in snp_annot_cols:
		#pearsonr( , )[0] is the statistic, [1] is the p-value
		cell_value = stats.pearsonr(gene_cols[cluster], snp_annot_cols[chunk])[0]
		#pval = stats.pearsonr(gene_cols[cluster], snp_annot_cols[chunk])[1]
		all_correlations.append(cell_value)

		if cluster == "021655":
			print(chunk)
			sys.exit()
		chunk1, chunk2 = chunk.split("-")
		phage2gene_1 = ""
		phage2gene_2 = ""
		for i in range(len(gene_coords)):
			if int(chunk1) >= gene_coords[i][0] and int(chunk1) <= gene_coords[i][1]:
				phage2gene_1 = i
			if int(chunk2) >= gene_coords[i][0] and int(chunk2) <= gene_coords[i][1]:
				phage2gene_2 = i
		if phage2gene_1 == "" and phage2gene_2 == "":
			phage2gene = chunk
		elif phage2gene_1 == phage2gene_2:
			phage2gene = "EF12" + str(76+phage2gene_1)
		else:
			if phage2gene_1 != "":
				a1 = "EF12" + str(76+phage2gene_1)
			if phage2gene_2 != "":
				a2 = "EF12" + str(76+phage2gene_2)
			phage2gene = a1 + "-" + a2
		
		epa_gene = ""
		if cluster in epa_orthogroups:
			epa_gene = epa_orthogroups[cluster]

		if cluster == "014989" or cluster == "022911":
			print("Correlation between", phage2gene, "and", epa_gene, "is", cell_value)

		if cluster in gene_functions:
			gene_function = gene_functions[cluster]
		else:
			gene_functions_local_list = []
			synerclust_file2 = open("results400/results/final_clusters.txt", "r")
			cluster_name = "Cluster" + ("0" * (6-len(cluster))) + cluster
			line = synerclust_file2.readline()
			while line != "":
				if line.startswith(cluster_name):
					genes = line.split("\t")[1].split()
					break
				line = synerclust_file2.readline()
			for gene in genes:
				genome = gene[0:-11]
				found = False
				for filename in os.listdir("annot"):
					if genome in filename:
						annot_file = open("annot/" + filename, "r")
						annot_line = annot_file.readline()
						while "##FASTA" not in annot_line:
							if "Prodigal" in annot_line and gene in annot_line:
								words = annot_line.split(";")
								for word in words:
									if "product=" in word:
										function_for_one_genome = word.replace("product=", "")
										gene_functions_local_list.append(function_for_one_genome)
										found = True
								break
							annot_line = annot_file.readline()
						break		
				if not found:
					print(gene, "was not found")
			most_common_gene_function = max(set(gene_functions_local_list), key = gene_functions_local_list.count)
			#If the number of times this gene is found is greater than half the length of the list. Checks if this function is found in the majority of the genomes
			if gene_functions_local_list.count(most_common_gene_function) > len(gene_functions_local_list)/2:
				gene_function = most_common_gene_function
				gene_functions[cluster] = gene_function
			elif (cell_value <= -0.6 or cell_value >= 0.6) and not printed:
				print(cluster_name + ":")
				for gene in set(gene_functions_local_list):
					print(gene + " count = " + str(gene_functions_local_list.count(gene)) + "\n")
				printed = True
		
		#if cell_value <= -0.6 or cell_value >= 0.6:
		#	output3.write(cluster + "\t" + gene_function.replace("\t", " ") + "\t" + phage2gene + "\t" + str(cell_value) + "\t" + epa_gene + "\n")

		if (int(chunk1) >= 11344 and int(chunk1) <= 13134) or (int(chunk2) >= 11344 and int(chunk2) <= 13134):
			EF1291_correlations.append(cell_value)
			#if cell_value <= -0.6 or cell_value >= 0.6:
			#	output4.write(cluster + "\t" + gene_function.replace("\t", " ") + "\t" + phage2gene + "\t" + str(cell_value) + "\t" + epa_gene + "\n")

#print("Max correlation =", max(saved))
#print("Min correlation =", min(saved))
"""histplot(all_correlations, binwidth=0.01)
plt.savefig("all_genes_correlation_histogram")
plt.show()
histplot(EF1291_correlations, binwidth=0.01)
plt.savefig("all_genes_with_1291_correlation_histogram")"""

#THIS IS THE MOST RECENT CODE
"""fig, ax = plt.subplots()
ax.hist(all_correlations, bins=100, log=True)
fig.savefig("all_genes_correlation_log_histogram")
fig2, ax2 = plt.subplots()
ax2.hist(EF1291_correlations, bins=100, log=True)
fig2.savefig("all_genes_with_1291_correlation_log_histogram")"""

#for tuple0 in sorted(saved, key=lambda item: item[0]):
#	print("Correlation between", tuple0[1], "and phage2 position", tuple0[2], "is", tuple0[0])
#sys.exit()
#End of columnize4"""

#Obsolete
#row_names = []
#for name in row_names2:
#	if int(name) > 12896 and int(name) < 13077:
#		row_names.append(name)

#df = DataFrame(data=corr_cols, index=row_names)
#clustermap(data=df, vmin=-1, vmax=1).savefig("orthogroup_correlation_heatmap")
#clustermap(data=df, row_cluster=False, vmin=-1, vmax=1).savefig("phage2_in_order_orthogroup_correlation_heatmap")
#sns.clustermap(data=df, vmin=-1, vmax=1).savefig("orthogroup_correlation_heatmap")
#sns.clustermap(data=df, row_cluster=False, vmin=-1, vmax=1).savefig("Phage2_in_order_orthogroup_correlation_heatmap")

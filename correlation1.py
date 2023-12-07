import sys
from scipy import stats
#import scipy
from pandas import DataFrame
import seaborn as sns
import matplotlib.pyplot as plt

#Dictionary of all column vectors made, cps, epa, and tag. Keys are the gene names (ex. EF_002456T0).
annot_cols = {}
#Having an arbitrary list of all the cps, epa, and tag genes and serotypes so that all columns stick to the same order
#row_names = []
def columnize(filename):
	file0 = open(filename, "r")
	line = file0.readline()
	while line != "DATA\n":
		line = file0.readline()
		if line.startswith("FIELD_LABELS"):
			genes_in_order = line.split()[1:]
	for i in range(len(genes_in_order)):
		if "cps" in filename:
			new_gene_name = "cps_" + genes_in_order[i]
		elif "tag" in filename:
			new_gene_name = "tag_" + genes_in_order[i]
		elif "epa" in filename:
			new_gene_name = "epa_" + genes_in_order[i]
		if genes_in_order[i] == "EF_002454T0":
			new_gene_name = "cpsF"
		genes_in_order[i] = new_gene_name
		annot_cols[new_gene_name] = []
	line = file0.readline()
	while line != "":
		for i in range(len(genes_in_order)):
			gene = genes_in_order[i]
			digit = int(line.split()[i+1])
			annot_cols[gene].append(digit)
		line = file0.readline()

#columnize("sorted_cps_heatmap.txt")
#columnize("sorted_epa_heatmap.txt")
#columnize("sorted_tag_heatmap.txt")
columnize("sorted_cps_heatmap400.txt")
columnize("sorted_epa_heatmap400.txt")
columnize("sorted_tag_heatmap400.txt")

#Uncomment to include capsule type columns
"""annot_cols["T1"] = []
annot_cols["T2"] = []
annot_cols["T5"] = []
serotype_file = open("sorted_serotype_colorstrip.txt", "r")
line = serotype_file.readline()
while line != "DATA\n":
	line = serotype_file.readline()
line = serotype_file.readline()
while line != "":
	serotype = line.split()[2]
	if serotype == "T1":
		annot_cols["T1"].append(1)
		annot_cols["T2"].append(0)
		annot_cols["T5"].append(0)
	elif serotype == "T2":
		annot_cols["T1"].append(0)
		annot_cols["T2"].append(1)
		annot_cols["T5"].append(0)
	elif serotype == "T5":
		annot_cols["T1"].append(0)
		annot_cols["T2"].append(0)
		annot_cols["T5"].append(1)
	else:
		annot_cols["T1"].append(0)
		annot_cols["T2"].append(0)
		annot_cols["T5"].append(0)
	line = serotype_file.readline()"""

snp_annot_cols = {}
def columnize2(filename):
	for i in range(1276, 1294):
		snp_annot_cols["EF" + str(i)] = []
	file0 = open("../nonsynonymous/" + filename, "r")
	if filename == "sorted_frameshift_only_heatmap.txt" or filename == "sorted_nonsyn_only_heatmap.txt":
		line = file0.readline()
		while line != "DATA\n":
			line = file0.readline()
	line = file0.readline()
	while line != "":
		for i in range(18):
			digit = int(line.split()[i+1])
			if digit == 0:
				snp_annot_cols["EF" + str(1276+i)].append(0)
			else:
				snp_annot_cols["EF" + str(1276+i)].append(1)
		line = file0.readline()

#Don't forget to uncomment the necessary lines after line 187
#columnize2("sorted_frameshift_only_heatmap.txt")
#columnize2("sorted_nonsyn_only_heatmap.txt")
#columnize2("sorted_nonsyn_and_no_frameshift_heatmap.txt")
#columnize2("sorted_nonsyn_and_frameshift_heatmap.txt")

def columnize3(filename):
	for i in range(1276, 1294):
		snp_annot_cols["EF" + str(i)] = []
	snp_annot_cols["EF1277-1278_intergenic"] = []
	snp_annot_cols["EF1283-1284_intergenic"] = []
	snp_annot_cols["EF1292-1293_intergenic"] = []
	file0 = open("../snps/" + filename, "r")
	line = file0.readline()
	while line != "DATA\n":
		line = file0.readline()
	line = file0.readline()
	while line != "":
		for i in range(18):
			digit = int(line.split()[i+1])
			if digit == 0:
				snp_annot_cols["EF" + str(1276+i)].append(0)
			elif digit == 1:
				snp_annot_cols["EF" + str(1276+i)].append(1)
			else:
				print("Not binary")
				sys.exit()
		digit = int(line.split()[19])
		if digit == 0:
			snp_annot_cols["EF1277-1278_intergenic"].append(0)
		elif digit == 1:
			snp_annot_cols["EF1277-1278_intergenic"].append(1)
		else:
			print("Not binary")
			sys.exit()
		
		digit = int(line.split()[20])
		if digit == 0:
			snp_annot_cols["EF1283-1284_intergenic"].append(0)
		elif digit == 1:
			snp_annot_cols["EF1283-1284_intergenic"].append(1)
		else:
			print("Not binary")
			sys.exit()
		
		digit = int(line.split()[21])
		if digit == 0:
			snp_annot_cols["EF1292-1293_intergenic"].append(0)
		elif digit == 1:
			snp_annot_cols["EF1292-1293_intergenic"].append(1)
		else:
			print("Not binary")
			sys.exit()
		line = file0.readline()

#Don't forget to uncomment the necessary lines after line 187
#columnize3("sorted_one_col_per_gene_binary_heatmap.txt")

def columnize4(filename):
	file0 = open(filename, "r")
	line = file0.readline()
	while line != "DATA\n":
		line = file0.readline()
		if line.startswith("FIELD_LABELS"):
			genes_in_order = line.split()[1:]
	for gene in genes_in_order:
		snp_annot_cols[gene] = []
	line = file0.readline()
	while line != "":
		for i in range(len(genes_in_order)):
			gene = genes_in_order[i]
			digit = int(line.split()[i+1])
			snp_annot_cols[gene].append(digit)
		line = file0.readline()
	return genes_in_order

#Don't forget to uncomment the necessary lines after line 187
#row_names = columnize4("../snps/sorted_150bp_binary_heatmap.txt")
row_names = columnize4("../snps/sorted_150bp_400genomes_heatmap.txt")
#row_names2 = columnize4("../snps/sorted_EF1291_heatmap_new.txt")
#row_names = columnize4("../nonsynonymous/sorted_nonsyn_EF1291_heatmap.txt")
#row_names = columnize4("../nonsynonymous/sorted_nonsyn_and_no_frameshift_EF1291_heatmap.txt")
#row_names = columnize4("../nonsynonymous/sorted_frameshift_only_EF1291_heatmap.txt")
#row_names = columnize4("../nonsynonymous/sorted_nonsyn_and_frameshift_EF1291_heatmap.txt")

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

"""If using columnize2:
corr_cols = {}
for capsule_gene in annot_cols:
	column = []
	for i in range(1276, 1294):
		if "EF" + str(i) not in phage2_genes_to_skip:
			#pearsonr( , )[0] is the statistic, [1] is the p-value
			cell_value = stats.pearsonr(annot_cols[capsule_gene], snp_annot_cols["EF" + str(i)])[0]
			column.append(cell_value)
	corr_cols[capsule_gene] = column

row_names = []
for i in range(1276, 1294):
	if "EF" + str(i) not in phage2_genes_to_skip:
		row_names.append("EF" + str(i))
#End of columnize2"""

"""If using columnize3:
corr_cols = {}
for capsule_gene in annot_cols:
	column = []
	for i in range(1276, 1294):
		if "EF" + str(i) not in phage2_genes_to_skip:
			#pearsonr( , )[0] is the statistic, [1] is the p-value
			cell_value = stats.pearsonr(annot_cols[capsule_gene], snp_annot_cols["EF" + str(i)])[0]
			column.append(cell_value)
	if "EF1277-1278_intergenic" not in phage2_genes_to_skip:
		cell_value = stats.pearsonr(annot_cols[capsule_gene], snp_annot_cols["EF1277-1278_intergenic"])[0]
		column.append(cell_value)
	if "EF1283-1284_intergenic" not in phage2_genes_to_skip:
		cell_value = stats.pearsonr(annot_cols[capsule_gene], snp_annot_cols["EF1283-1284_intergenic"])[0]
		column.append(cell_value)
	if "EF1292-1293_intergenic" not in phage2_genes_to_skip:
		cell_value = stats.pearsonr(annot_cols[capsule_gene], snp_annot_cols["EF1292-1293_intergenic"])[0]
		column.append(cell_value)	
	corr_cols[capsule_gene] = column

row_names = []
for i in range(1276, 1294):
	if "EF" + str(i) not in phage2_genes_to_skip:
		row_names.append("EF" + str(i))
if "EF1277-1278_intergenic" not in phage2_genes_to_skip:
	row_names.append("EF1277-1278_intergenic")
if "EF1283-1284_intergenic" not in phage2_genes_to_skip:
	row_names.append("EF1283-1284_intergenic")
if "EF1292-1293_intergenic" not in phage2_genes_to_skip:
	row_names.append("EF1292-1293_intergenic")
#End of columnize3"""

"""If using columnize4:"""
#for finding the positions of the blackest cells in EF1291
#saved = []
saved = {}

corr_cols = {}
for capsule_gene in annot_cols:
	column = []
	for chunk in snp_annot_cols:
		#if int(chunk) > 12896 and int(chunk) < 13077:
		#pearsonr( , )[0] is the statistic, [1] is the p-value
		cell_value = stats.pearsonr(annot_cols[capsule_gene], snp_annot_cols[chunk])[0]
		pval = stats.pearsonr(annot_cols[capsule_gene], snp_annot_cols[chunk])[1]
		column.append(cell_value)
		if abs(cell_value) > 0.5:
			print("Phage 2 =", chunk, "cell wall gene =", capsule_gene, "correlation =", cell_value, "p-value =", pval)
			#if capsule_gene in ["epa_EF_002160T0", "epa_EF_002161T0", "epa_EF_002162T0", "epa_EF_002163T0", "epa_EF_002164T0"]:
				#print("Phage 2 =", chunk, "epa =", capsule_gene, "correlation =", cell_value, "p-value =", pval)
			#else:
				#print("ALERT", "Phage 2 =", chunk, "epa =", capsule_gene, "correlation =", cell_value, "p-value =", pval)
#		if capsule_gene == "epa_EF_002164T0":
#			print(chunk, "correlation =", cell_value)
		#if cell_value < -0.4:
			#saved.append(cell_value)
			#saved[chunk] = cell_value
		#	saved[chunk] = (cell_value, pval)
	corr_cols[capsule_gene] = column
#for finding the positions of the blackest cells in EF1291
#sns.histplot(saved, binwidth=0.01)
#plt.savefig("correlation_distribution_nonsyn")

#for tuple0 in sorted(saved.items(), key=lambda item: item[1][0]):
	#print(tuple0[0] + "\t" + str(tuple0[1]))
#	print(tuple0[0] + "\t" + str(tuple0[1][0]) + "\t" + str(tuple0[1][1]))
#sys.exit()
#End of columnize4"""

#Only for the 3rd call to columnize4:
#row_names = []
#for name in row_names2:
#	if int(name) > 12896 and int(name) < 13077:
#		row_names.append(name)

#df = DataFrame(data=corr_cols, index=row_names)
#Only uncomment one of the following three
#sns.clustermap(data=df, vmin=-1, vmax=1).savefig("150bp_400_genomes_correlation_heatmap")
#sns.clustermap(data=df, row_cluster=False, vmin=-1, vmax=1).savefig("150bp_400_genomes_correlation_heatmap_in_order")
#sns.clustermap(data=df, metric="correlation").savefig("correlation_heatmap_diff_metric")

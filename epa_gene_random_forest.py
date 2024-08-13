import sys
from sklearn.ensemble import RandomForestClassifier as rf
import numpy as np
import copy

epa_gene = "000091"
number_of_random_forests = 100
topX = 50
percent = 0.8

"""Building the SNP input matrix"""
snp_dict = {}
genomes_file = open("../operon/all_epa_on_one_contig_genomes.txt", "r")
line = genomes_file.readline()
while line != "":
	genome = line.replace("\n", "")
	if "006494875" not in genome and "902164645" not in genome:
		snp_dict[line.replace("\n", "")] = []
	line = genomes_file.readline()

rep_genomes = []
sequences = []
alignment = open("/gsap/archive-bacterial/Projects/EnteroGenome/phage2/all_faecalis_whole_genome_phylogeny/snps/rep_alignment.afa", "r")
line = alignment.readline()
while line != "":
	if not line.startswith(">"):
		print("Error")
		sys.exit()
	genome = line[1:-1]
	sequence = alignment.readline().replace("\n", "")
	rep_genomes.append(genome)
	sequences.append(sequence)
	line = alignment.readline()
num_genomes = len(rep_genomes)
ref_index = rep_genomes.index("GCA_006494835.1_ASM649483v1")
ref_sequence = sequences[ref_index]
all_positions = []
for g in range(num_genomes):
	seq = sequences[g]
	row = []
	for i in range(len(ref_sequence)):
		if seq[i] != ref_sequence[i]:
			row.append(1)
		else:
			row.append(0)
	all_positions.append(row)
diff_columns = []
for i in range(len(all_positions[0])):
	all_same = True
	for g in range(num_genomes):
		if all_positions[g][i] == 1:
			all_same = False
			break
	if not all_same:
		diff_columns.append(i)
#print("Number of SNP positions in EF1291:", len(diff_columns))
num_complete_epa_found = 0
once = True
for g in range(num_genomes):
	genome = rep_genomes[g]
	if genome in snp_dict:
		num_complete_epa_found += 1
		for i in diff_columns:
			snp_dict[genome].append(all_positions[g][i])
		if once and genome != "GCA_006494835.1_ASM649483v1":
			#print(genome)
			#print("Make sure this isn't all 0's:", snp_dict[genome])
			once = False
#print("Number of rep alignment genomes found in complete epa file:", num_complete_epa_found)

X = []
order_of_genomes_used_in_ml = []
for genome in snp_dict:
	order_of_genomes_used_in_ml.append(genome)
	X.append(snp_dict[genome])

"""Building the epa presence labels (the "correct answers" for whether each genome has the epa gene or not)"""
orthogroups_in_order = []
y = []
#for i in range(len(order_of_genomes_used_in_ml)):
#	y.append([])
orthogroup_file = open("../operon/synerclust/results6/results/cluster_dist_per_genome.txt", "r")
header = orthogroup_file.readline()
line = orthogroup_file.readline()
while not line.startswith(epa_gene):
	line = orthogroup_file.readline()
#while line != "":
	#if "0" in line.split()[2:]:
		#orthogroups_in_order.append(line.split()[0])
		#i = 0
for genome in order_of_genomes_used_in_ml:
	col = header.split().index(genome + ".fa")
	y.append(line.split()[col])
			#y[i].append(line.split()[col])
	#i += 1
	#line = orthogroup_file.readline()
#print(y)

#clf1 = rf(oob_score = True, n_estimators=200)
#clf2 = rf(oob_score = True, n_estimators=200, max_samples=0.75)
#clf3 = rf(oob_score = True, n_estimators=200, max_samples=0.7)
#clf4 = rf(oob_score = True, n_estimators=200, max_samples=0.65)
#clf5 = rf(oob_score = True, n_estimators=200, max_samples=0.2)
#clf6 = rf(oob_score = True, max_features=40)
#clf7 = rf(oob_score = True, max_features=55)
#clf1.fit(X, y)
#clf2.fit(X, y)
#clf3.fit(X, y)
#clf4.fit(X, y)
#clf5.fit(X, y)
#clf6.fit(X, y)
#clf7.fit(X, y)
#print(clf1.oob_score_)
#print(clf2.oob_score_)
#print(clf3.oob_score_)
#print(clf4.oob_score_)
#print(clf5.oob_score_)
#print(clf6.oob_score_)
#print(clf7.oob_score_)
#print("Num features:", clf.n_features_in_)
#print("Classes:", clf1.classes_)
#print("n classes:", clf1.n_classes_)
#print("Outputs:", clf1.n_outputs_)
#print(clf1.feature_importances_)

#Run random forest once and print out the top 10 SNPs
"""results = copy.deepcopy(clf1.feature_importances_)
for i in range(10):
	index = np.argmax(results)
	print("Position", diff_columns[index], "in alignment has highest feature importance", max(results))
	results = np.delete(results, index)
	del diff_columns[index]
"""
#Do random forest multiple times
top50_counts = {}
average_importance = {}
first = True
for i in range(number_of_random_forests):
	clf1 = rf(oob_score = True, n_estimators=200)
	clf1.fit(X, y)
	counter = 0
	for i, val in sorted(enumerate(clf1.feature_importances_), key=lambda x:x[1], reverse=True):
		if counter < topX:
			if i in top50_counts:
				top50_counts[i] += 1
			else:
				top50_counts[i] = 1
		if i in average_importance:
			average_importance[i] += val
		else:
			average_importance[i] = val
		counter += 1
#print(top50_counts)
for i in top50_counts:
	if top50_counts[i] >= percent * number_of_random_forests:
		print(str(diff_columns[i]) + ": average feature importance=" + str(average_importance[i]/number_of_random_forests))

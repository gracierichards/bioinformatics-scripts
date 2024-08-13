import sys
from sklearn.cluster import AgglomerativeClustering as cluster
#Uncomment if you need the lines that print the version number
#import sklearn
#Imports for dendrogram graphing
import numpy as np
#from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram

#print(sklearn.__version__)
#sys.exit()

user_defined_threshold = 6

#List of genomes to use in itol in order
#itol_genomes_output = open("genomes_grouped_by_epa_8clusters.txt", "w")
colorstrip = open("threshold_" + str(user_defined_threshold) + "_colorstrip.txt", "w")
color_dict = {0:"#ff0000", 1:"#ff8000", 2:"#ffff00", 3:"#00ff00", 4:"#00b000", 5:"#00ffff", 6:"#0000ff", 7:"#a000ff", 8:"#00a0ff", 9:"#ffb0e0", 10:"#a0a0a0", 11:"#000000", 12:"#800000", 13:"#007000", 14:"#706000", 15:"#c080ff", 16:"#ffc0a0", 17:"#ffc000", 18:"#c0ff00", 19:"#00ffd0"}

data = open("results6/results/cluster_dist_per_genome.txt", "r")
header = data.readline()
genomes_in_order = header.split()[2:]
i = genomes_in_order.index("GCA_000007785.1_ASM778v1.fa.ref")
line = data.readline()
#X = np.array([[]] * (len(genomes_in_order)-1))
matrix = []
for j in range(len(genomes_in_order)-1):
	matrix.append([])
while line != "":
	row = 0
	digits = line.split()
	for j in range(2, len(digits)):
		#Deleting the ref genome column because it's not in the SNP file
		if j-2 != i:
			digit = digits[j]
			if digit == "0":
#				np.append(X[row], 0)
				matrix[row].append(0)
				row += 1
			elif digit == "1":
#				np.append(X[row], 1)
				matrix[row].append(1)
				row += 1
			elif digit == "2":
#				np.append(X[row], 1)
				matrix[row].append(1)
				row += 1
			else:
				print("Error: Cell is not 0, 1, or 2. Value:", digit)
	line = data.readline()

print("Number of rows is", len(matrix))
print("Number of columns is", len(matrix[0]))

genomes_in_order.remove("GCA_000007785.1_ASM778v1.fa.ref")

def cluster_func(threshold):
	print("distance_threshold=" + str(threshold))
	clustering = cluster(distance_threshold=threshold, n_clusters=None).fit(matrix)
	num_in_each_cluster = {}
	for val in clustering.labels_:
		if val in num_in_each_cluster:
			num_in_each_cluster[val] += 1
		else:
			num_in_each_cluster[val] = 1
	print(num_in_each_cluster)
	
	genomes_in_each_cluster = {}
	for i in range(len(genomes_in_order)):
		cluster_name = clustering.labels_[i]
		if cluster_name in genomes_in_each_cluster:
			genomes_in_each_cluster[cluster_name].append(genomes_in_order[i])
		else:
			genomes_in_each_cluster[cluster_name] = [genomes_in_order[i]]
		colorstrip.write(genomes_in_order[i].replace(".fa", "") + " " + color_dict[cluster_name] + " " + str(cluster_name) + "\n")

	for cluster_name in genomes_in_each_cluster:
		print("Last genome in cluster is", genomes_in_each_cluster[cluster_name][-1])
		#for genome in genomes_in_each_cluster[cluster_name]:
			#itol_genomes_output.write(genome + "\n")

	if len(clustering.labels_) != len(genomes_in_order):
		print("Lengths are different,", len(genomes_in_order), "genomes were included but the length of the results list is", len(clustering.labels_))

cluster_func(user_defined_threshold)

#clustering = cluster(n_clusters=10).fit(matrix)
"""clustering = cluster(distance_threshold=0, n_clusters=None).fit(matrix)
clusters = {}
for val in clustering.labels_:
	if val in clusters:
		clusters[val] += 1
	else:
		clusters[val] = 1
print(clusters)
if len(clustering.labels_) != len(genomes_in_order):
	print("Lengths are different,", len(genomes_in_order), "genomes were included but the length of the results list is", len(clustering.labels_))

def plot_dendrogram(model, **kwargs):
	counts = np.zeros(model.children_.shape[0])
	n_samples = len(model.labels_)
	for i, merge in enumerate(model.children_):
		current_count = 0
		for child_idx in merge:
			if child_idx < n_samples:
				current_count += 1  # leaf node
			else:
				current_count += counts[child_idx - n_samples]
		counts[i] = current_count
	linkage_matrix = np.column_stack( [model.children_, model.distances_, counts] ).astype(float)
	dendrogram(linkage_matrix, **kwargs)

plt.title("Dendrogram from Hierarchical Clustering Epa Types")
#plot_dendrogram(clustering, truncate_mode="level", p=3)
plot_dendrogram(clustering)
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.savefig("epa_dendrogram")"""

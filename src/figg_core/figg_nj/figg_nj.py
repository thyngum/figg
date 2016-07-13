import logging

"""
Functions to build the neighbor-joining tree
"""

def symetric_val(matrix, i, j):
	"Returns the symetric value of the matrix if operations move towards the upper triangle"

	if j > i:
		i, j = j, i
			
	return matrix[i][j]


def total_genome_dist(matrix, i):
	"Returns the sums of distances from one genome to all the others"

	dist = sum([symetric_val(matrix, i, j) for j in range(len(matrix))])

	return dist


def dist_norm(matrix, i, j):
	"Normalizes the distance between two genomes"
	
	# Takes the distance between two genomes and subtracts the total distances of these
	# two genomes from all other genomes scaled by (number of genomes - 2).
	dist = (symetric_val(matrix, i, j) - \
		   (1.0/(len(matrix) - 2))*(total_genome_dist(matrix, i) + total_genome_dist(matrix, j)))

	return dist


def edge_length(matrix , i, j):
	"Calculates distance between new combined genome and old genomes"

	dist = (0.5*(symetric_val(matrix, i, j) + \
	       (1.0/(len(matrix) - 2))*(total_genome_dist(matrix, i) - total_genome_dist(matrix, j))))

	return dist


def new_matrix_val(matrix, z, i, j):
	"Calculates the distances between the combine genome at i,j with the other genomes at k"

	dist = 0.5*(symetric_val(matrix, i, z) + symetric_val(matrix, j, z) - symetric_val(matrix, i, j))

	return dist


def nj(matrix, genomes, heights):

	if ( len(matrix) == 2 ):
		# Base case or last iteration of the NJ algorithm

		# Tree heights
		heights += [matrix[1][0]/2]

		# The actual NJ tree in Newick format
		tree =  "(" + genomes[0] + ", " + genomes[1] + ": " + str(matrix[1][0]) + ");"

		return tree, heights

	else:
		minimum = 0
		index_i = 0
		index_j = 0
		for i in range(len(matrix)):
			temp = min([dist_norm(matrix, i, j) for j in range(len(matrix)) if i != j])
			if min(temp ,minimum) == temp:
				index_j = [dist_norm(matrix , i, j) for j in range(len(matrix)) if i != j].index(temp)
				index_i = i
				minimum = temp
			
		tree = "(" + str(genomes[index_i]) + ": " + str(edge_length(matrix,index_i,index_j)) +\
			", " + str(genomes[index_j]) + ": " + str(edge_length(matrix,index_j,index_i)) + ")"		
		
		new_genomes = genomes[:]
		new_genomes.remove(genomes[index_i])
		new_genomes.remove(genomes[index_j])
		new_genomes.insert(0,tree)
		heights += [matrix[index_i][index_j]/2]

		new_matrix = [[0]*(len(matrix) - 1) for x in range(len(matrix)-1)]
		for z in range(1, len(new_matrix)):
			ind = genomes.index(new_genomes[z])
			new_matrix[z][0] = new_matrix_val(matrix, ind, index_i, index_j)

		for k in range(2, len(new_matrix)):
			for g in range(1, k):
			    ind_i = genomes.index(new_genomes[k])
			    ind_j = genomes.index(new_genomes[g])
			    new_matrix[k][g] = matrix[ind_i][ind_j]

		return nj(new_matrix, new_genomes, heights)

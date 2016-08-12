#!/usr/bin/env python

import argparse
import logging
import os
import sys

# Import FIGG core modules
import figg_core.figg_parser as figg_parser			# Functions to parse input format
import figg_core.figg_matrices as figg_matrices		# Functions to compute distance matrices
import figg_core.figg_output as figg_output			# Functions to write output

# import upgma.upgma as upgma 						# Basic implementation of the UPGMA algorithm (TODO)
import nj.nj as nj									# Basic implementation of the NJ algorithm

def main():
    
	# Configure logging
	logging.basicConfig(level = logging.INFO)
        
	# Define input parameters
	parser = argparse.ArgumentParser(description = 'figg')
	parser.add_argument('input_file', help = 'input file')
	parser.add_argument('--circular', dest = 'circular', action = 'store_const',
	                    default = False, const = True, help = 'flag indicating that input genomes are circular')
	parser.add_argument('--output-format', dest = 'output_format', default = 'text',
	 					help = 'output format, can be either \'text\' (default), \'mega\', \'nexus\', \'phylip\' or \'all\'') 
	parser.add_argument('--verbose', dest = 'verbose', action = 'store_const',
	                    default = False, const = True,  help = 'show the results of all runtime calculations')
	args = parser.parse_args()

	# Read input parameters
	input_file = args.input_file
	is_circular = args.circular
	output_format = args.output_format
	supported_formats = ['text', 'mega', 'nexus', 'phylip', 'all' ]
	if ( output_format not in supported_formats ):
		sys.exit("Format '%s' not supported (valid formats are 'text', 'mega', 'nexus', 'phylip' or 'all')" % output_format)
	verbose = args.verbose
	if ( verbose ):
	    if ( is_circular ):
	        genome_type = 'circular'
	    else:
	        genome_type = 'linear'    

	    # logging.info("Input file: %s (genomes are %s)" % ( input_file, genome_type ))   
	    print "Input file: %s (genomes are %s)" % ( input_file, genome_type )  
	    # logging.info("Output format: %s" % output_format) 
	    print "Output format: %s" % output_format   

	# Initialize global variables
	genome_labels = []                # List of genome labels
	gene_orders = []                  # List with the order of genes in each input genome
	num_genomes = 0                   # Number of genomes in input
	ref_order = []		              # Reference gene order
	ref_matrix = []                   # Reference adjacency matrix
	first_seen = []                   # Genome in which each gene was first seen
	adj_matrices = []                 # Set of adjacency matrices for each genome in input
	dist_matrix = []                  # Matrix of observed differences
	corrected_dist_matrix = []        # Matrix of corrected differences

	# Read input file
	genome_labels, gene_orders = figg_parser.parse_file(input_file, is_circular)
	num_genomes = len(gene_orders)

	# Define the workspace (referece order and matrix)
	ref_order, ref_matrix, first_seen = figg_matrices.workspace(gene_orders, genome_labels, is_circular)

	# Compute adjacency matrices and distances
	adj_matrices = figg_matrices.adj_matrix_set(ref_matrix, ref_order, gene_orders)
	dist_matrix = figg_matrices.dist_matrix(adj_matrices)
	pos_freq_matrix = figg_matrices.freq_matrix_pos(adj_matrices)
	neg_freq_matrix = figg_matrices.freq_matrix_neg(adj_matrices)
	corrected_dist_matrix = figg_matrices.dist_matrix_corrected(adj_matrices, pos_freq_matrix, neg_freq_matrix)

	'''
	# Build the NJ tree 
	# tree, heights = nj.tree(dist_matrix, genome_labels, [])
	# tree_corrected, heights_corrected = nj.tree(corrected_dist_matrix, genome_labels, [])
	'''

 	# Write the output
	abs_path = os.path.abspath(input_file)
	output_path, input_filename = os.path.split(abs_path)
	output_prefix = output_path + '/' + '.'.join(input_filename.split('.')[:-1])

	if ( output_format == 'text' or output_format == 'all' ):
		figg_output.print_matrix_to_file(dist_matrix, output_prefix + "_dist_matrix.tsv", genome_labels)	
		figg_output.print_matrix_to_file(corrected_dist_matrix, output_prefix + "_corrected_dist_matrix.tsv", genome_labels)

	if ( output_format == 'mega' or output_format == 'all' ):
		figg_output.print_mega_format(dist_matrix, genome_labels, output_prefix + "_dist_matrix.meg")	
		figg_output.print_mega_format(corrected_dist_matrix, genome_labels, output_prefix + "_corrected_dist_matrix.meg")

	if ( output_format == 'phylip' or output_format == 'all' ):
		figg_output.print_phylip_format(dist_matrix, genome_labels, output_prefix + "_dist_matrix.phy")	
		figg_output.print_phylip_format(corrected_dist_matrix, genome_labels, output_prefix + "_corrected_dist_matrix.phy")	

	if ( output_format == 'nexus' or output_format == 'all' ):
		figg_output.print_nexus_format(dist_matrix, genome_labels, output_prefix + "_dist_matrix.nex")	
		figg_output.print_nexus_format(corrected_dist_matrix, genome_labels, output_prefix + "_corrected_dist_matrix.nex")

	if ( verbose ):

		# Print basic info
		print "Number of genomes: %i" % num_genomes
		print "Number of genes in the workspace: %i" % len(ref_order)
		print "  Gene\tFirst seen in"
		for i in range(len(ref_order)):
			print "  %s\t%s" % (ref_order[i], first_seen[i])
		print ""

		# Print adjacency matrices for each genome
		print "Adjacency matrices for each genome:\n"
		for i in range(num_genomes):
			figg_output.print_matrix(adj_matrices[i], ref_order)		

		# Print frequency matrices
		print "Positive and negative frequency matrices:\n"
		figg_output.print_matrix(pos_freq_matrix, ref_order)
		figg_output.print_matrix(neg_freq_matrix, ref_order)

		# Print distance matrices
		print "Uncorrected distance matrix:\n"		
		figg_output.print_matrix(dist_matrix, genome_labels)
		print "Corrected distance matrix:\n"
		figg_output.print_matrix(corrected_dist_matrix, genome_labels)

		'''
		# Print the NJ tree
		print "NJ tree from uncorrected distance matrix:\n"
		print tree, "\n"
		print "NJ tree from corrected distance matrix:\n"
		print tree_corrected, "\n"
		'''

		print "Output saved to '%s'." % output_path


if __name__ == '__main__':
	try: 
		main()
	except Exception, error:
		logging.error("%s" % str(error))	
    
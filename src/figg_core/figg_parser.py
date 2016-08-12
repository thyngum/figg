"""
Functions to parse input format
"""

def parse_file(input_file, is_circular):

	labels = []
	orders = []

	# Read input file
	f = open(input_file, "rU")
	temp = f.read()
	temp = temp.split('>')
	del temp[0]

	# Parse the input
	for i in range(len(temp)):
		labels.append(temp[i].split('\n')[0])
		orders.append(temp[i].split('\n')[1])
	for i in range(len(orders)):
		orders[i] = orders[i].split(" ")

	# If genomes are circular then append the first gene at the end of each gene list
	if is_circular:
		[i.append(i[0]) for i in orders]

	return labels, orders

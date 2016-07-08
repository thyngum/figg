#!/usr/bin/env python

import argparse
import sys
import logging

import figg_core.figg_core as figg_core

def main():
    
    # configures logging
    logging.basicConfig(level=logging.INFO)
        
    # Parses input parameters
    parser = argparse.ArgumentParser(description='figg')
    parser.add_argument('input_file', help='Input file')
    parser.add_argument('--circular', dest='circular', action='store_const',
                        default=False, const=True, 
                        help='A flag indicating that input genomes are circular.')
    parser.add_argument('--output-format', dest='output_format', action='store', 
                        required = True, choices=set(["mega", "paup", "phylip", "all"]),
                        help='Indicate the output format (mega, paup, phylip, all)')
    parser.add_argument('--verbose', dest='verbose', action='store_const',
                        default=False, const=True, 
                        help='Output the results of all runtime calculations.')
    args = parser.parse_args()
    
    # Read input parameters
    input_file = args.input_file
    is_circular = args.circular
    output_format = args.output_format
    verbose = args.verbose

    # logging.info("Input file: [%s]" % input_file)
    # logging.info("Circular genomes: [%s]" % str(is_circular))
    # logging.info("Output format: [%s]" % output_format)
    if ( verbose ):
        if ( is_circular ):
            genome_type = 'circular'
        else:
            genome_type = 'linear'    
        print "Input file: %s (genomes are %s)" % ( input_file, genome_type )   
        print "Output format: %s" % output_format   
    
    # Calls the program
    figg_core.run_figg(input_file, is_circular, output_format, verbose)

"""    
    try:
        figg_core.run_figg(input_file, is_circular, output_format)
        logging.info("Done!")
    except Exception, e:
        logging.error("Figg raised an error: [%s]" % str(e))
"""    
    
if __name__ == '__main__':
    main()
    
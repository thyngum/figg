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
                        help='A flag indicating if the genome is circular.')
    parser.add_argument('--output-format', dest='output_format', action='store', 
                        required = True, choices=set(["mega", "paup", "phylip", "all"]),
                        help='The output format (mega, paup, phylip, all)')
    args = parser.parse_args()
    
    # Read input parameters
    input_file = args.input_file
    logging.info("Input file: [%s]" % input_file)
    is_circular = args.circular
    logging.info("Circular genomes: [%s]" % str(is_circular))
    output_format = args.output_format
    logging.info("Output format: [%s]" % output_format)
    
    # Calls the program
    figg_core.run_figg(input_file, is_circular, output_format)

"""    
    try:
        figg_core.run_figg(input_file, is_circular, output_format)
        logging.info("Done!")
    except Exception, e:
        logging.error("Figg raised an error: [%s]" % str(e))
"""    
    
if __name__ == '__main__':
    main()
    
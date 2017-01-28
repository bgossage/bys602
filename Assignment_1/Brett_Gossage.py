#!/usr/bin/python
#
#
# BYS 602, Spring 2017
#
# Brett Gossage (bgossgen@gmail.com)
#
# Assignment 1, Jan 27, 2017
#


import sys

#
# Setup path to genebank module...
#
sys.path.append( "../modules" )

import genebank

filename = "S_invicta_sequence.gb"

# If provided, use the filename from the command line arg...
if( len(sys.argv) > 1 ):
   filename = sys.argv[1] 
   
# Create a GenBank file parser...
parser = genebank.GenBankParser()

# Define the set of keywords we are going to process...
filter_keys = { "DEFINITION", "ACCESSION", "AUTHORS", "CDS", "gene" }

parser.set_filtered( filter_keys )

# Create a list of gene locus data records to be filled by the parser...
loci = []

# Open the output file and assure is gets closed...
with open( "seq_parse.txt", "w" ) as output_file:

# Parse the given GenBank file...
   parser.parse( filename, loci )

# Print out the list of locus records...
   for locus in loci:
   # Print the locus to the screen...
      print locus
   # Also write the request data to the output file...
      output_file.write( "\n\nLOCUS: \n" )
      
      output_file.write( "Definition: {:s}\n".format( locus.definition ) )
      output_file.write( "Accession: {:s}\n".format( locus.accession ) )
      output_file.write( "Authors: {:s}\n".format( locus.authors ) )
      output_file.write( "CDS start and stop locations: {:d} to {:d}\n".format(locus.cds_location[0],locus.cds_location[1]) )
      output_file.write( "Gene start and stop locations: {:d} to {:d}\n".format(locus.cds_location[0],locus.cds_location[1]) )
      
      output_file.write( "///\n" )
      
#end with output_file


# EOF

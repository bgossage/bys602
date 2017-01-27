#!/usr/bin/python
#
#
# BYS 602, Spring 2017
#
# Brett Gossage (bgossgen@gmail.com)
#
# Assignment 1, Jan 14, 2017
#


import sys

#
# Setup path to genebank module...
#
sys.path.append( "../modules" )

import genebank

filename = "S_invicta_sequence.gb"


if( len(sys.argv) > 1 ):
   filename = sys.argv[1] 
   
parser = genebank.GenBankParser()

filter_keys = { "DEFINITION", "ACCESSION", "AUTHORS", "CDS", "gene" }

parser.set_filtered( filter_keys )

loci = []

with open( "seq_parse.txt", "w" ) as output_file:

   parser.parse( filename, loci )

   for locus in loci:
      
      print locus
      
      output_file.write( "\n\nLOCUS: \n" )
      
      output_file.write( "Definition: {:s}\n".format( locus.definition ) )
      output_file.write( "Accession: {:s}\n".format( locus.accession ) )
      output_file.write( "Authors: {:s}\n".format( locus.authors ) )
      output_file.write( "CDS start and stop locations: {:d} to {:d}\n".format(locus.cds_location[0],locus.cds_location[1]) )
      output_file.write( "Gene start and stop locations: {:d} to {:d}\n".format(locus.cds_location[0],locus.cds_location[1]) )
      
      output_file.write( "///\n" )
      
#end with output_file

   


# EOF

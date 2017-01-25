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

filter_keys = { "DEFINITION" }

parser.set_filtered( filter_keys )

with open( "seq_parse.txt", "w" ) as output_file:

     parser.parse( filename, output_file )
     
#end with output_file
   


# EOF

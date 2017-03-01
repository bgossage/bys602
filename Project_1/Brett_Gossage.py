#!/usr/bin/python
#
#
# BYS 602, Spring 2017
#
# Brett Gossage (bgossgen@gmail.com)
#
# Project 1, March 3, 2017
#


import sys
import termcolor

#
# Setup path to genalign module...
#
sys.path.append( "../modules" )

import genalign
import genebank


# Get the filenames from the command line...
if( len(sys.argv) < 3 ):
    print "Usage: ", sys.argv[0], " <genebank file> <sequence file>"
    sys.exit( -1 )
    
genebank_filename = sys.argv[1]
sequence_filename = sys.argv[2]

gap_penalty = -2.0
  
similarityMatrix = genalign.SimilarityMatrix()
 
similarityMatrix.readFrom( "JN_blosom50.txt" )

# Create a GenBank file parser...
parser = genebank.GenBankParser()

# Create a list of gene locus data records to be filled by the parser...
loci = []

# Parse the given GenBank file...
parser.parse( genebank_filename, loci )

# Read the sequence fasta file..
seq = genebank.ReadFasta( sequence_filename )

for locus in loci:
    
    scoreMatrix = genalign.ScoringMatrix( similarityMatrix, seq, locus.origin, gap_penalty )
      
    aligned_seq1, aligned_seq2 = scoreMatrix.backtrace( seq, loci.origin )

    print "Score: ", scoreMatrix.score()
    print "Aligned Seq1: ", termcolor.colored( aligned_seq1, "red" )
    print "Aligned Seq2: ", termcolor.colored( aligned_seq1, "blue" )
    print "--------------------------------------\n"
    
#end for locus
    

   



# EOF

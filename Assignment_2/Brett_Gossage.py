#!/usr/bin/python
#
#
# BYS 602, Spring 2017
#
# Brett Gossage (bgossgen@gmail.com)
#
# Assignment 2, Feb 11, 2017
#


import sys

#
# Setup path to genalign module...
#
sys.path.append( "../modules" )

import genalign

seq1 = "IQIFSFIFRQEWNDA"  
seq2 = "QIFFFFRMSVEWND" 

print "Seq1: ", seq1
print "Seq2: ", seq2

match = 5.0
mismatch = -1.0
gap_penalties = [ 0, -2, -4, -8, -16 ]
  
similarityMatrix = genalign.SimilarityMatrix(match,mismatch)  

def RunCases( matrix_file ):
 
    similarityMatrix.readFrom( "JN_blosom50.txt" )
    
    for gap_penalty in gap_penalties:
      
      scoreMatrix = genalign.ScoringMatrix( similarityMatrix, seq1, seq2, gap_penalty )
      
      aligned_seq1, aligned_seq2 = scoreMatrix.backtrace( seq1, seq2 )

      print "Similarity Matrix: ", matrix_file
      print "Gap penalty: ", gap_penalty

      print "Aligned Seq1: ", aligned_seq1
      print "Aligned Seq2: ", aligned_seq2
      print "--------------------------------------\n"
   

RunCases( "JN_blosom50.txt" )

RunCases( "BNG_blosom80.txt" )

RunCases( "BNG_PAM10.txt" )

RunCases( "BNG_PAM80.txt" )

# EOF

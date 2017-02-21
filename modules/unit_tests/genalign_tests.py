#!/usr/bin/python

"""
GeneBank module unit tests.

Each class/method is tested against expected results from the given inputs.

(Uses the unit test module from the python libraries (similar to JUnit.)

"""

import sys

sys.path.append( ".." )


import unittest
import genalign


class GeneAlign_tests( unittest.TestCase ):

   """ 
      .
   """
   def test_sequence_align( self ):
              
      seq1 = "GCATGCU"  
      seq2 = "GATTACA"   
      
      match = 1.0
      mismatch = -1.0
      gap_penalty = -1.0
              
      similarityMatrix = genalign.SimilarityMatrix(match,mismatch,"GCTAU")
      
      scoreMatrix = genalign.ScoringMatrix( similarityMatrix, seq1, seq2, gap_penalty )

      print scoreMatrix.score_matrix
      print "score = ", scoreMatrix.score()
      #print scoreMatrix.arrow
           
      scoreMatrix.backtrace( seq1, seq2 )
      
      aligned_seq1, aligned_seq2 = scoreMatrix.backtrace( seq1, seq2 )
      
      print "Seq1: ", seq1
      print "Seq2: ", seq2
      print "Aligned Seq1: ", aligned_seq1
      print "Aligned Seq2: ", aligned_seq2
      
      self.assertEqual( aligned_seq1, "GCATG_CU" )
      self.assertEqual( aligned_seq2, "G_ATTACA" )
      
   #/////////////////////////////   
      seq1 = "IQIFSFIFRQEWNDA"  
      seq2 = "QIFFFFRMSVEWND" 
      
      similarityMatrix.readFrom( "JN_blosom50.txt" )
      
      scoreMatrix = genalign.ScoringMatrix( similarityMatrix, seq1, seq2, gap_penalty )
      
      #print scoreMatrix.score_matrix
      #print scoreMatrix.arrow
      
      aligned_seq1, aligned_seq2 = scoreMatrix.backtrace( seq1, seq2 )

      print "Seq1: ", seq1
      print "Seq2: ", seq2
      print "Aligned Seq1: ", aligned_seq1
      print "Aligned Seq2: ", aligned_seq2

   #end test_sequence_align ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   """ 
      .
   """
   def test_read_matrix( self ):
       
       similarityMatrix = genalign.SimilarityMatrix()
       
       similarityMatrix.readFrom( "JN_blosom50.txt" )
       
       self.assertEqual( similarityMatrix.indexMap[0], "A" )
       self.assertEqual( similarityMatrix.indexMap[22], "Z" )
       
       self.assertEqual( similarityMatrix.matrix[0,0], 5.0 )
       self.assertEqual( similarityMatrix.matrix[22,22], 5.0 )
       self.assertEqual( similarityMatrix.matrix[23,23], 1.0 )
       
   #end test_read_matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
   def test_write_matrix( self ):
       
       match = 1.0
       mismatch = -1.0

       similarityMatrix = genalign.SimilarityMatrix(match,mismatch,"GCTAU")
       
       similarityMatrix.writeTo( "simil.txt" )
       
   #end test_write_matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# end class GeneBank_tests

if __name__ == '__main__':
    unittest.main()

# EOF

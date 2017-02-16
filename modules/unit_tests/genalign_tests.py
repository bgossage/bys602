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
              
      seq1 = "CGCA"  
      seq2 = "CACGCAT"   
      
      match = 1.0
      not_match = 0.0
      gap_penalty = -1.0
              
      substMatrix = genalign.SubstitutionMatrix(match,not_match)
      
      scoreMatrix = genalign.ScoringMatrix( substMatrix, seq1, seq2, gap_penalty )

      #print scoreMatrix.score_matrix
      #print scoreMatrix.arrow
           
      scoreMatrix.backtrace( seq1, seq2 )
      
      aligned_seq1, aligned_seq2 = scoreMatrix.backtrace( seq1, seq2 )
      
      print seq1
      print seq2
      print aligned_seq1
      print aligned_seq2
      
   #/////////////////////////////   
      seq1 = "IQIFSFIFRQEWNDA"  
      seq2 = "QIFFFFRMSVEWND" 
      
      substMatrix.readFrom( "JN_blosom50.txt" )
      
      scoreMatrix = genalign.ScoringMatrix( substMatrix, seq1, seq2, gap_penalty )
      
      #print scoreMatrix.score_matrix
      #print scoreMatrix.arrow
      
      aligned_seq1, aligned_seq2 = scoreMatrix.backtrace( seq1, seq2 )

      print seq1
      print seq2
      print aligned_seq1
      print aligned_seq2

     # self.assertEqual( loci[0].definition, "Solenopsis invicta venom allergen 3 (LOC105199703), mRNA." )
      
   #end test_sequence_align ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   """ 
      .
   """
   def test_read_matrix( self ):
       
       substMatrix = genalign.SubstitutionMatrix()
       
       substMatrix.readFrom( "JN_blosom50.txt" )
       
       self.assertEqual( substMatrix.indexMap[0], "A" )
       self.assertEqual( substMatrix.indexMap[22], "Z" )
       
       self.assertEqual( substMatrix.matrix[0,0], 5.0 )
       self.assertEqual( substMatrix.matrix[22,22], 5.0 )
       self.assertEqual( substMatrix.matrix[23,23], 1.0 )
       
   #end test_read_matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# end class GeneBank_tests

if __name__ == '__main__':
    unittest.main()

# EOF

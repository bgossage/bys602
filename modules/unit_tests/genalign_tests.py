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
import numpy
import math

class GeneAlign_tests( unittest.TestCase ):

   """ 
      .
   """
   def test_sequence_align( self ):
              
      seq1 = "ARND"  
      seq2 = "ARNE"   
      
      match = 1.0
      not_match = -1.0
              
      substMatrix = genalign.SubstitutionMatrix(match,not_match)
      
      scoreMatrix = genalign.ScoringMatrix( substMatrix, seq1, seq2 )

      print scoreMatrix.score_matrix
      print scoreMatrix.arrow
           
      scoreMatrix.backtrace( seq1, seq2 )
      
      aligned_seq1, aligned_seq2 = scoreMatrix.backtrace( seq1, seq2 )
      
      print aligned_seq1
      print aligned_seq2

     # self.assertEqual( loci[0].definition, "Solenopsis invicta venom allergen 3 (LOC105199703), mRNA." )
      
   #end test_parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# end class GeneBank_tests

if __name__ == '__main__':
    unittest.main()

# EOF

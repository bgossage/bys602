#!/usr/bin/python

"""
GeneBank module unit tests.

Each class/method is tested against expected results from the given inputs.

(Uses the unit test module from the python libraries (similar to JUnit.)

"""

import sys

sys.path.append( ".." )


import unittest
import genebank
import numpy
import math

class GeneBank_tests( unittest.TestCase ):

   """ 
      Test parsing a GenBank file.
   """
   def test_parse( self ):
      
      parser = genebank.GenBankParser()
      
      filter_keys = { "DEFINITION", "ACCESSION" }

      parser.set_filtered( filter_keys )

      loci = []
      
      parser.parse( "invicta_small.gb", loci )


      self.assertEqual( loci[0].definition, "Solenopsis invicta venom allergen 3 (LOC105199703), mRNA." )
      self.assertEqual( loci[0].accession, "NM_001304591 XM_011166905" )

   #end test_parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# end class GeneBank_tests

if __name__ == '__main__':
    unittest.main()

# EOF

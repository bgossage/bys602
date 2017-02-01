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

      loci = []
      
      parser.parse( "invicta_small.gb", loci )

      for locus in loci:
         print locus

      self.assertEqual( loci[0].definition, "Solenopsis invicta venom allergen 3 (LOC105199703), mRNA." )
      self.assertEqual( loci[0].accession, "NM_001304591 XM_011166905" )
      self.assertEqual( loci[0].authors, "Schmidt M, McConnell TJ and Hoffman DR." )
      self.assertEqual( loci[0].cds_location[0], 1 )
      self.assertEqual( loci[0].cds_location[1], 705 )
      
      self.assertEqual( loci[1].definition, "Solenopsis invicta queen venom protein Sol i IV precursor, mRNA,\n\
            partial cds." )
      self.assertEqual( loci[1].accession, "AY963564" )
      self.assertEqual( loci[1].authors, "Deslippe,R.J., HaghiPour-Peasley,J.D., San Francisco,M.J., Fokar,M.\n\
            and Renthal,R.D." )
      self.assertEqual( loci[1].cds_location[0], 1 )
      self.assertEqual( loci[1].cds_location[1], 408 )
      self.assertEqual( loci[1].origin[0], 'a' )
      
   #end test_parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# end class GeneBank_tests

if __name__ == '__main__':
    unittest.main()

# EOF

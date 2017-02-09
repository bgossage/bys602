#!/usr/bin/python
#
#
# BYS 602, Spring 2017
#
# Brett Gossage (bgossgen@gmail.com)
#
# 
#

"""@package genalign

A module for computing sequence alignments.

"""

import sys
import re
import string
import numpy

##
# 
#
class  SubstitutionMatrix:

##
## Constructor.
#
   def __init__( self ):

      self.matrix = numpy.identity(20,float)
      
      self.indexMap = "ARNDCQEGHILKMFPSTWYV"

   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~
   
## Conversion to string for printing.
##
   def __str__( self ):
  
      ans = str()
      
      return ans
      
   #end str ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#end class SubstitutionMatrix

##
# .
#
class  ScoringMatrix:

##
## Constructor.
#
   def __init__( self, substMatrix, seq1, seq2, gapPenalty=2.0 ):
      
      length1 = len(seq1) + 1
      length2 = len(seq2) + 1
      
      shape = length1, length2
   
      self.score_matrix = numpy.zeros( shape, numpy.int )
      
      self.arrow = numpy.zeros( shape, numpy.int )
      
   # Fill the first row and first column...
      self.score_matrix[0] = numpy.arange(length2) * gapPenalty
      self.score_matrix[:,0] = numpy.arange(length1) * gapPenalty
      
   # Fill the first row of the arrow array...
      self.arrow[0] = numpy.ones(length2)

   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~
   
## Conversion to string for printing.
##
   def __str__( self ):
  
      ans = str()
      
      return ans
      
   #end str ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#end class ScoringMatrix


## EOF

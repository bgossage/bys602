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

   def compare( self, aa1, aa2 ):

      n1 = self.indexMap.index( aa1 )
      n2 = self.indexMap.index( aa2 )

      return self.matrix[n1,n2]

   #endif compare ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
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
      
      f = numpy.zeros(3)

   # 
      for i in range(1,length1):
         for j in range(1,length2):
            
            f[0] = self.score_matrix[i-1,j] + gapPenalty
            f[1] = self.score_matrix[i,j-1] + gapPenalty
            f[2] = self.score_matrix[i,j-1] + substMatrix.compare(seq1[i-1],seq2[j-1])

            self.score_matrix[i,j] = f.max()

            self.arrow[i,j] = f.argmax()
            
         #end for j
      #end for i
      
   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~
   
   def backtrace( self, seq1, seq2 ):
   
      str1, str2 = "", "" # The resulting alignment as a pair of strings
      
      v,h = self.arrow.shape
      
      ok = True
      v -= 1
      h -= 1
      
   #   while ok:
         
         
      #end while
   
   #end backtrace() ~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
## Conversion to string for printing.
##
   def __str__( self ):

      ans = str()

      return ans

   #end str ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#end class ScoringMatrix


## EOF

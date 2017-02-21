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
import string
import StringIO

##
# A Similarity matrix class
#
class  SimilarityMatrix:

##
## Constructor.
# @param match the match score
# @param not_match the mismatch score
#
   def __init__( self, match=1.0, mismatch=0.0, indexMap="ARNDCQEGHILKMFPSTWYV" ):

      size = len(indexMap)

      self.matrix = numpy.ndarray( (size,size),float)
      
      self.matrix.fill( mismatch )
      
      numpy.fill_diagonal( self.matrix, match )
      
      self.indexMap = indexMap

   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~

   ##
   # Read a substitution matrix from a file.
   # @param filename the name of the input file.
   #
   def readFrom( self, filename ):
   #
   # Open the input file and cleanly handle exceptions...
   #
      with open( filename, "r" ) as file:
         
         header_read = False
         size = 0
         row = 0         
         
         for line in file:

         # Skip comments...
             if line.startswith("#"): continue
                 
             line = line.strip()
             tokens = line.split()
             
         # Skip empty lines...
             if not tokens: continue
                 
         # Look for the index map row...
             if not header_read:
                 self.indexMap = string.join( tokens, "" )
                 size = len(self.indexMap)
                 self.matrix = numpy.ndarray( (size,size ),float )
                 header_read = True
             else:
             # Store current row data from the line values...
                 for i in range(0,size):
                     self.matrix[row,i] = float(tokens[i+1])
                 row += 1
         #end for   
          
      #end with
    
   #end readFrom ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ##
   # Lookup the substitution likelihood parameter for a given 
   # pair of sequence elements
   #
   def similarity( self, aa1, aa2 ):

      n1 = self.indexMap.index( aa1 )
      n2 = self.indexMap.index( aa2 )

      return self.matrix[n1,n2]

   #endif compare ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
## Convert a SimilarityMatrix to string for printing.
##
   def writeTo( self, filename ):
   
      matrix_str = numpy.array2string( self.matrix )
      
      with open( filename, "w" ) as file:
         file.write( "# A similarity matrix\n" )
         file.write( self.indexMap +"\n" )
         file.write( matrix_str )


      
   #end str ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#end class SimilarityMatrix

##
# A scoring matrix using the Needleman-Wunseh method.
#
class  ScoringMatrix:
# A sequence alignment scoring matrix.
##
## Constructor.
#
   def __init__( self, similarityMatrix, seq1, seq2, gapPenalty=-2.0 ):
   ##
   # Create a scoring matrix using the Needleman-Wunsch method.
   #
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
             
            prob = similarityMatrix.similarity(seq1[i-1],seq2[j-1])  
              
            f[0] = self.score_matrix[i-1,j] + gapPenalty
            f[1] = self.score_matrix[i,j-1] + gapPenalty
            f[2] = self.score_matrix[i-1,j-1] + prob

            self.score_matrix[i,j] = f.max()

            self.arrow[i,j] = f.argmax()
            
         #end for j
      #end for i
      
   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~
   

   def backtrace( self, seq1, seq2 ):
   ## Backtrace this scoring matrix to align the given sequences.
      str1, str2 = "", "" # The resulting alignment as a pair of strings
      
      v,h = self.arrow.shape
      
      ok = 1
      v -= 1
      h -= 1
       
      while ok:
         direction = self.arrow[v,h]

         if direction == 0: # left
            str1 += seq1[v-1]
            str2 += "_"
            v -= 1
         elif direction == 1: # up
            str1 += "_"
            str2 += seq2[h-1]
            h -= 1
         elif direction == 2: # diagonal
            str1 += seq1[v-1]
            str2 += seq2[h-1]
            v -= 1
            h -= 1
         if v <= 0 or h <= 0:
             ok = 0
             
      # end while
             
   #reverse the strings...
      str1 = str1[::-1]
      str2 = str2[::-1]
      
      return str1, str2
   
   #end ScoringMatrix.backtrace() ~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   def score( self ):
   
      rows, cols = self.score_matrix.shape
      
      return self.score_matrix[rows-1,cols-1]
      
   #end ScoringMatrix.score() ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
## Conversion to string for printing.
##
   def __str__( self ):

      ans = str()

      return ans

   #end str ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#end class ScoringMatrix


## EOF

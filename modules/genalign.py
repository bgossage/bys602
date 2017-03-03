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
import re
import termcolor

##
# A sequence, either nucleic or protein
#
class Sequence:
    
    def __init__(self, name="", sequence=""):
        self.name= name
        self.data = sequence 
    #end __init__
        
    def length( self ):
        return len(self.data)
        
    ## Overload operator []
    def __getitem__( self, index ):
        return self.data[index]
        
    ## Overload operator [] 
    def __setitem__( self, index, value ):
        self.data[index] = value  
        
#end class Sequence

##
# Read a fasta file and return a list of FASTA records.
#
def ReadFasta( filename ):
    
   with open( filename, "r" ) as fasta_file:

       sequence_list = []  # the returned sequence list
       sequence = None
       record_num = -1
       
       for line in fasta_file:
           
           if not line: continue
           
           if line.startswith(">"):
               sequence = Sequence( name=line[1:])
               sequence_list.append( sequence )
               record_num += 1
           else:
               if record_num > -1:
               # We are reading a sequence.  Add subset to the sequence...
                  subseq = re.sub(r'\s+', '', line)  # Remove whitespace
                  sequence.data += subseq
           
       #end for
       
       return sequence_list
       
   #end with   
    
#end ReadFasta ~~~~~~~~~~~~~~~~~~~~~~~~~~

##
# The alignment of two sequences.
#
class Alignment:
    
    def __init__( self, query_seq=Sequence(), subj_seq=Sequence() ):
        
        self.name = ""
        
        self.query_seq = query_seq
        self.subj_seq = subj_seq
        
        self.score = 0.0
  
    #end Alignment.__init__  
    
    def write( self, outfile ):
        
        outfile.write( "Query: {:s}\n".format(self.query_seq.name ))
        outfile.write( "Subj: {:s}\n".format(self.subj_seq.name ))
        outfile.write( "score: {:6.2f}\n\n".format(self.score))
        
        block_length = 60        
        count = 0
        length = self.subj_seq.length()
        
        while count < length:
            
            if count + block_length > length:
                block_length = length % block_length
            
            outfile.write( "Subj:  {:4d}  {:s}\n".format( count+1,
                                                          self.subj_seq[count:count + block_length])
                                                        )
            
            outfile.write( "Query: {:4d}  ".format( count+1 ) ) 
        
            for i in range(count,count + block_length):
                char = self.query_seq[i]
                if char != self.subj_seq[i]:
                   outfile.write( termcolor.colored( char, "red" ) )
                else:
                   outfile.write(char)
            #end for
            outfile.write( "\n\n" )
            
            count += block_length
                   
        #end while          
        
    #end write ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
#end class Alignment
        

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
      try:
         n1 = self.indexMap.index( aa1 )
         n2 = self.indexMap.index( aa2 )
      except(ValueError) as error:
         msg = str(error) + " for keys: " + aa1 + ":" + aa2
         raise Exception(msg)

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
   def __init__( self, similarityMatrix, query_seq, subj_seq, gapPenalty=-2.0 ):
   ##
   # Create a scoring matrix using the Needleman-Wunsch method.
   #
      length1 = query_seq.length() + 1
      length2 = subj_seq.length() + 1
      
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
             
            prob = similarityMatrix.similarity(query_seq[i-1],subj_seq[j-1])  
              
            f[0] = self.score_matrix[i-1,j] + gapPenalty
            f[1] = self.score_matrix[i,j-1] + gapPenalty
            f[2] = self.score_matrix[i-1,j-1] + prob

            self.score_matrix[i,j] = f.max()

            self.arrow[i,j] = f.argmax()
            
         #end for j
      #end for i
      
   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~
   

   def backtrace( self, query_seq, subj_seq ):
   ## Backtrace this scoring matrix to align the given sequences.
   
      ans = Alignment() # The resulting alignments
      
      v,h = self.arrow.shape
      
      ok = 1
      v -= 1
      h -= 1
       
      while ok:
         direction = self.arrow[v,h]

         if direction == 0: # left
            ans.query_seq.data += query_seq[v-1]
            ans.subj_seq.data += "_"
            v -= 1
         elif direction == 1: # up
            ans.query_seq.data += "_"
            ans.subj_seq.data += subj_seq[h-1]
            h -= 1
         elif direction == 2: # diagonal
            ans.query_seq.data += query_seq[v-1]
            ans.subj_seq.data += subj_seq[h-1]
            v -= 1
            h -= 1
         if v <= 0 or h <= 0:
             ok = 0
      # end while
             
   #reverse the strings...
      ans.query_seq = Sequence( name=query_seq.name, sequence=ans.query_seq[::-1] )
      ans.subj_seq = Sequence( name=subj_seq.name,   sequence=ans.subj_seq[::-1] )
      
      ans.score = self.score()
      
      return ans
   
   #end ScoringMatrix.backtrace() ~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   def score( self ):
   
      rows, cols = self.score_matrix.shape
      
      return self.score_matrix[rows-1,cols-1]
      
   #end ScoringMatrix.score() ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
## Conversion to string for printing.
##
   def __str__( self ):

      ans = str()
      
      ans = numpy.array2string( self.score_matrix )
      
      return ans

   #end str ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#end class ScoringMatrix


## EOF

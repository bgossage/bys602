#!/usr/bin/python
#
#
# BYS 602, Spring 2017
#
# Brett Gossage (bgossgen@gmail.com)
#
# 
#

import sys

class  GenBankLocus:

#
# Constructor.
#
   def __init__( self ):
      
      self.definition = str()
      
      self.accession = str()

      self.authors = str()

   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~
   
   def __str__( self ):
      
      return "Locus\n Definition: {:s}\n Accession: {:s}\n".format(self.definition, self.accession)
      
   #end str ~~~~~~~~~~~~~~~~~~~~~~~~

#end class GenBankLocus


def handle_definition( record, locus ):

   loc = record.find( "DEFINITION" )

   data = record[loc+12 : len(record)-1]

   locus.definition = data

#end handle_definition() ~~~~~~~~~~~~~~~~~~~~~~~~

def handle_accession( record, locus ):

   loc = record.find( "ACCESSION" )

   data = record[loc+12 : len(record)-1]

   locus.accession = data

#end handle_definition() ~~~~~~~~~~~~~~~~~~~~~~~~


class GenBankParser:
   
#
# Constructor.
#
   def __init__( self ):

      self.keywords = {
                         "DEFINITION", "AUTHORS", "ACCESSION",
                         "TITLE", "REFERENCE", "FEATURES", "ORIGIN",
                         "LOCUS", "JOURNAL", "ORGANISM", "VERSION",
                         "KEYWORDS", "source", "gene", "CDS", 
                         "sig_peptide", "mat_peptide" 
                         "//"
                      }

      self.handlers = dict()

      self.handlers["DEFINITION"] = handle_definition
      self.handlers["ACCESSION"] = handle_accession

      self.filtered = { }

   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_filtered( self, filter_set ):

      self.filtered = filter_set

   #end set_filtered ~~~~~~~~~~~~~~~~~~~~

   def parse( self, filename, loci ):

      reading_record = False
      key = ""
      read_key = ""
      content = ""
      locus = 0

      #
      # Open the input file and cleanly handle exceptions.
      #
      with open( filename, "r" ) as genbank_file:

         for line in genbank_file:

            if line.startswith( "LOCUS" ):
               locus = GenBankLocus()
               loci.append(locus)

            tokens = line.split()

            if not tokens:
               continue

            key = tokens[0]

         # If we are reading a record for a given key...
            if reading_record :
            # And we encounter any keyword..
               if key in self.keywords:
               # We've reached the end of the record...
                  reading_record = False

               # Process the record content...
                  if read_key in self.handlers:
                     self.handlers[read_key]( content, locus )

               # Clear the content string...
                  content = ""
               #end if key
            # end if reading

            # Otherwise, add the current line to the record content...
               else:
                  content += line

            #end if

         # If we are not reading a record...
            if not reading_record:

            # And if we encounter a key to be read...
               if key in self.filtered:
               # We've reached the beginning of a record... 
                  reading_record = True
                  read_key = key

               # Add the current line to the record content...
                  content += line
            #end if

         # If we encounter the end of a locus,
            if key == "//":

               loci.append( locus )

            #end if

         # end for line

      # end with genbank_file

   #end parse ~~~~~~~~~~~~~~~~~~~~~~

#end class GenBankParser


class GenBankRecordHandler:
   
#
# Constructor.
#
   def __init__( self, line ):
      
      tokens = line.split()

      self.key = tokens[0]

      self.content = line[ len(key):len(line)-1 ]

   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~
   
   
   def handle( self, line ):
   
      self.content += line
 
   # end handle
   
   
   
#end class GenBankRecordHandler

## EOF

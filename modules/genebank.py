#!/usr/bin/python
#
#
# BYS 602, Spring 2017
#
# Brett Gossage (bgossgen@gmail.com)
#
# 
#

"""@package genebank
A module for parsing and storing GenBank data.

"""

import sys
import re

##
# A GenBank locus record
#
class  GenBankLocus:

##
# Constructor.
#
   def __init__( self ):

      self.definition = str()
      
      self.accession = str()

      self.authors = str()
      
      self.cds_location = 0,0
      
      self.gene_location = 0,0

   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~
   
## Conversion to string for printing.
#
   def __str__( self ):
  
      ans = "Locus\n Definition: {:s}\n".format(self.definition, self.accession)
      ans += " Accession: {:s}\n".format(self.accession)
      ans += " Authors: {:s}\n".format(self.authors)
      ans += " CDS location: {:d}..{:d}\n".format(self.cds_location[0],self.cds_location[1])
      ans += " Gene location: {:d}..{:d}\n".format(self.gene_location[0],self.gene_location[1])
      
      return ans
      
   #end str ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#end class GenBankLocus

## Extract content from a GenBank Record
#
# @param key the record key e.g. AUTHORS
# @param record the record content lines
#
def extract_content( key, record ):
   
   loc = record.find( key )

   data = record[loc+len(key) : ]

   return data.strip()

#end extract_content ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def extract_range( data ):

   if (data.find( "join" ) >= 0) or (data.find( "complement" ) >= 0) :
      pair = 0,0
      return pair

   end = data.find("\n")

   begin = 0
   if data[0] == "<":
      begin += 1

   loc = data[begin:end]

   elip = loc.find("..")
   
   if( elip < 0 ):
      pair = 0,0
      return pair

   start = int(loc[0:elip])

   if loc[elip+2] == ">":
      elip += 1

   stop =  int(loc[elip+2:])

   pair = start, stop
   
   return pair

#end extract_range ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def handle_definition( record, locus ):

   data = extract_content( "DEFINITION", record )

   locus.definition = data

#end handle_definition() ~~~~~~~~~~~~~~~~~~~~~~~~

def handle_accession( record, locus ):

   data = extract_content( "ACCESSION", record )

   locus.accession = data

#end handle_accession() ~~~~~~~~~~~~~~~~~~~~~~~~

def handle_authors( record, locus ):

   data = extract_content( "AUTHORS", record )

   locus.authors = data

#end handle_authors() ~~~~~~~~~~~~~~~~~~~~~~~~

def handle_cds( record, locus ):

   data = extract_content( "CDS", record )

   locus.cds_location = extract_range( data )

#end handle_cds() ~~~~~~~~~~~~~~~~~~~~~~~~

def handle_gene( record, locus ):

   data = extract_content( "gene", record )

   locus.gene_location = extract_range( data )

#end handle_gene() ~~~~~~~~~~~~~~~~~~~~~~~~

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
      self.handlers["AUTHORS"] = handle_authors
      self.handlers["CDS"] = handle_cds
      self.handlers["gene"] = handle_gene

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
      ignore = False # Ignore text where keywords should not be read e.g. comments

   #
   # Open the input file and cleanly handle exceptions.
   #
      with open( filename, "r" ) as genbank_file:

         for line in genbank_file:

            if line.startswith( "LOCUS" ):
               locus = GenBankLocus()
               
         # Ignore words that occur within comments...
            if line.count( "\"" ) % 2 == 1:
                ignore = not ignore

            tokens = line.split()

            if not tokens:
               continue

            key = tokens[0]


         # If we are reading a record for a given key...
            if reading_record :
            # And we encounter any keyword..
               if key in self.keywords and not ignore:
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
            # Add the current locus to the list...
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

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
import string

##
# A GenBank locus record
#
class  GenBankLocus:

##
## Constructor.
#
   def __init__( self ):

      self.definition = str()
      
      self.accession = str()

      self.authors = str()
      
      self.cds_location = 0,0
      
      self.gene_location = 0,0
      
      self.origin = str()

   # end constructor ~~~~~~~~~~~~~~~~~~~~~~~~
   
## Conversion to string for printing.
##
   def __str__( self ):
  
      ans = "Locus\n Definition: {:s}\n".format(self.definition)
      ans += " Accession: {:s}\n".format(self.accession)
      ans += " Authors: {:s}\n".format(self.authors)
      ans += " CDS location: {:d}..{:d}\n".format(self.cds_location[0],self.cds_location[1])
      ans += " Gene location: {:d}..{:d}\n".format(self.gene_location[0],self.gene_location[1])
      ans += " Origin: {:s}\n".format(self.origin)
      
      return ans
      
   #end str ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#end class GenBankLocus

## Extract content for a given key from a GenBank record
##
## @param key the record key e.g. AUTHORS
## @param next_key the following key
## @param data the locus record content lines
##
## @return a string containing the content for the key
##
def extract_content( key, next_key, data ):

   regex = ".*"        # Skip everything up to the key
   regex += key        # Match the key
   regex += r"\s*"     # Skip white space
   regex += r"(.*?)"   # Capture group 1: characters up to the next key and don't be greedy
   regex += next_key

   pattern = re.compile( regex,
                         re.DOTALL # Include line feeds.
                       )

   loc = pattern.match( data )

   if None == loc:
      print "No Key: ", key
      return ""

   content = loc.group(1).strip()

   content.rstrip("\n")

   return content

#end extract_content ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Extract a range string for a given key from a GenBank record
##
## @param key the rance key e.g. CDS
## @param data the locus record content lines
##
## @return a pair of integers defining the range
##
def extract_range( key, data ):

   regex = ".*"       # Skip everything up to the key
   regex += key       # Match the key
   regex += r"\s*"    # Skip white space
   regex += r"<?"     # Match optional '<'
   regex += r"(\d+)"  # Capture group 1: one more digits
   regex += r"..>?"   # Elipsis with optional '>'
   regex += r"(\d+)"  # Capture group 2: one more digits

   pattern = re.compile( regex,
                         re.DOTALL # Include line feeds.
                       )

   loc = pattern.match( data )
   
   if None == loc:
#      print "No range for key: ", key
      pair = 0,0
      return pair
   
   start = int(loc.group(1))
   stop =  int(loc.group(2))

   pair = start, stop
   
   return pair
   
#end extract_range ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def parse_definition( record, locus ):

   data = extract_content( key="DEFINITION", next_key="ACCESSION", data=record )

   locus.definition = data

#end parse_definition() ~~~~~~~~~~~~~~~~~~~~~~~~

def parse_accession( record, locus ):

   data = extract_content( key="ACCESSION", next_key="VERSION", data=record )

   locus.accession = data

#end parse_accession() ~~~~~~~~~~~~~~~~~~~~~~~~

def parse_authors( record, locus ):

   data = extract_content( key="AUTHORS", next_key="TITLE", data=record )

   locus.authors = data

#end parse_authors() ~~~~~~~~~~~~~~~~~~~~~~~~

def parse_cds( record, locus ):

   locus.cds_location = extract_range( key="CDS", data=record )

#end parse_cds() ~~~~~~~~~~~~~~~~~~~~~~~~

def parse_gene( record, locus ):

   locus.gene_location = extract_range( key="gene", data=record )

#end parse_gene() ~~~~~~~~~~~~~~~~~~~~~~~~

def parse_origin( record, locus ):
    
# Find the end of the ORIGIN key...
   origin_pattern = re.compile( r"\bORIGIN\b" )
   origin_match = origin_pattern.search( record )
   start = origin_match.span()[1]  # start here to find the origin data

# The sequence is broken into blocks separated white space and may be nucleotides or amino acids
# 
   nucleic_pattern = re.compile( r"(\b[acgturykmswbdhvn-]+\b)", re.IGNORECASE )
   
   amino_pattern = re.compile( r"(\b[ABCDEFHHIKLMNPQRSTUVWYZX*-]+\b)", re.IGNORECASE )

# Search for a nucleic sequence...
   sequence = nucleic_pattern.findall( record, start )

# If that fails, search for a protein sequence...
   if None == sequence:
      seqs = amino_pattern.findall( record, start )
    
# If both fail just return...
   if None == sequence:
       return
       
# Join the blocked data into a single string...
   locus.origin = string.join( sequence, "" )

#end parse_origin() ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Parse a record and store it in a locus object
##
def parse_locus( record, locus ):

   parse_origin( record, locus )

   parse_definition( record, locus )

   parse_accession( record, locus )

   parse_authors( record, locus )

   parse_cds( record, locus )

   parse_gene( record, locus )

#end parse_locus() ~~~~~~~~~~~~~~~~~~~~~~~~

##
## A parser for GenBack records (loci)
##
class GenBankParser:

## Parse and store a GenBank data file.
##
   def parse( self, filename, loci ):

      content = ""
      locus = 0

   #
   # Open the input file and cleanly handle exceptions.
   #
      with open( filename, "r" ) as genbank_file:
         
         for line in genbank_file:

            if line.startswith( "LOCUS" ):
               locus = GenBankLocus()

            content += line

         # If we encounter the end of a locus,
            if line.startswith("//") :
            # Parse the locus data...
               parse_locus( content, locus )

            # Add the current locus to the list...
               loci.append( locus )

            # Clear the content string...
               content = ""

            #end if

         # end for line

      # end with genbank_file

   #end parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#end class GenBankParser

def ReadFasta( filename ):
    
   with open( filename, "r" ) as fasta_file:
       
       lines = []
       
       for line in fasta_file:
           
           if line.startswith(">"): 
               continue

       # Replace all whitespace...
           data = re.sub(r'\s+', '', line)
               
           lines.append( data )
           
       #end for
           
       ans = string.join( lines, "" )
       
       return ans
   #end with   
    
#end ReadFasta ~~~~~~~~~~~~~~~~~~~~~~~~~~

## EOF

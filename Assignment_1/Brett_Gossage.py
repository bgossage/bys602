#!/usr/bin/python
#
#
# BYS 602, Spring 2017
#
# Brett Gossage (bgossgen@gmail.com)
#
# Assignment 1, Jan 14, 2017
#


import sys


filename = "S_invicta_sequence.gb"

keywords = {
              "DEFINITION", "AUTHORS", "ACCESSION",
              "TITLE", "REFERENCE", "FEATURES", "ORIGIN",
              "LOCUS", "JOURNAL", "ORGANISM", "VERSION",
              "KEYWORDS", "source", "CDS", "//"
           }

filtered = { "DEFINITION", "AUTHORS", "ACCESSION", "CDS" }

reject = keywords - filtered

print reject

if( len(sys.argv) > 1 ):
   filename = sys.argv[1] 
   
reading_record = False
key = ""
content = ""

#
# Open the input file and cleanly handle exceptions.
#
with open( filename, "r" ) as seq_file:

#
# Open the output file and cleanly handle exceptions., "\n"
#
   with open( "seq_parse.txt", "w" ) as output_file:

      line_num = 0
      
      for line in seq_file:
         
         if line.startswith( "LOCUS" ):
            print "======== START RECORD ============================"
            output_file.write( "======== START RECORD =============================\n" )

         tokens = line.split()

         if not tokens:
            continue
         
         key = tokens[0]
         
      # If we are reading a record for a given key...
         if reading_record :
         # And we encounter a keyword, we've reached the end of the record...
            if key in keywords:
            # Transition the parser state...
               reading_record = False

            # Process the record content...
               print content
               output_file.write( content )

            # Clear the content string...
               content = ""
         # end if
         
         # Otherwise, add the current line to the record content...
            else:
               content += line

         if not reading_record:
         # If we are reading a record...
            if key in filtered:
            # Transition the parser state...
               reading_record = True
            # Add the current line to the record content...
               content += line
         #end if

         if tokens[0] == "//":
            print "======== END RECORD =============================\n"
            output_file.write( "======== END RECORD =============================\n\n" )
         #end if

      # end for line

   # end with outputfile

# end with seq_file

# EOF

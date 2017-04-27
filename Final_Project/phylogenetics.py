#! /usr/local/python

import os
from Bio import SeqIO
from Bio import Phylo
import matplotlib
import pylab

try:
    
    input_filename = "S_invicta_mtDNA_query.nwk"

    tree = None

    with open(input_filename, "r" ) as Infile:
        
        tree = Phylo.read( Infile, "newick" )
    #end with

    tree.ladderize()
    
    pylab.cla()  #clear any prior graph    
    
    Phylo.draw( tree, branch_labels=lambda c: c.branch_length )
    
    pylab.xlabel( "clade" )
    pylab.show()
    
    with open("tree.txt","w") as Outfile:
        
       Phylo.draw_ascii( tree, Outfile )
       
    #end with
    
except Exception as error:
    print "Error!"
    print str(error)
    
   
# EOF



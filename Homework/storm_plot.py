#! /usr/local/python

from Bio import SeqIO
#import matplotlib
import pylab

Infile = open("Pol_d_venom_duped.fasta")
record_iterator = SeqIO.parse(Infile, "fasta")
rec_one = next(record_iterator)
rec_two = next(record_iterator)
Infile.close

seq_one = str(rec_one.seq).upper()
seq_two = str(rec_two.seq).upper()


window = 5 # Look for sub-sequences of this length
dict_one = {}
dict_two = {}

print "Seq1: ", seq_one
print "Seq2: ", seq_two

for (seq, section_dict) in [ (seq_one, dict_one),(seq_two, dict_two)]:
   for i in range(len(seq)- window):
        section = seq[i:i+window] 
        # If the section has already been searched
        if section_dict.has_key( section ):
        # Append the current index...
           section_dict[section].append(i) # At this 
        else:
        # Found a new sub-sequence. Added a new key/value pair.
        # NOTE: The section is the key and the value is a list of indices.
           section_dict[section] = [i] 

# Do not use exceptions for flow control
#	   try: #try-except block to handle exceptions
#			section_dict[section].append(i)
#	   except KeyError:
#			section_dict[section] = [i]
			
# Find any sub-sequences found in both sequences
# Note: dictionary iteration is over the keys only, so
#       this is the intersection of the subseqences.
matches = set(dict_one).intersection(dict_two)
print("%i unique matches" %len(matches) )

#Create lists of x and y co-ordinates for scatter plot

x = []
y = []
for section in matches:
	for i in dict_one[section]:
		for j in dict_two[section]:
			x.append(i)
			y.append(j)
			
pylab.cla()  #clear any prior graph
pylab.gray()
pylab.scatter(x,y)
pylab.xlim(0, len(seq_one) - window)
pylab.ylim(0, len(seq_two) - window)
pylab.xlabel("%s (length %i aa)" % (rec_one.id, len(seq_one) ) )
pylab.ylabel("%s (length %i aa)" % (rec_two.id, len(seq_two) ) )
pylab.show()
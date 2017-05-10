
HYPOTHESIS:

Eusocial behavior in ants and wasps developed first in a common ancestor. Molecular analysis of genomes or proteins from ants and eusocial wasps will show a higher degree of relatedness than that between ants and asocial wasp genera.

NOTE:  Current theory holds that nest building preceded eusocial behavior and all three genera examined
       for this project are nest builders.
       
PROCEDURE:

Sequence Acquisition:

To limit the amount of data to search, representative genera and species from each group were selected
as follows:

   Solitary Wasp: Trypoxylon politum (Organ pipe mud dauber)

   Social Wasp: Polistes dominulus (European paper wasp)

   Ant: Solenopsis invicta: (Red imported fire ant)
   
Mitochondrial DNA (mtDNA) partial coding sequences for cytochrome oxidase subunit I (COI)
for each species was obtained via NCBI search.
   
A protein BLAST search was done through NCBI using a S.invicta COI protein sequence as the query.
The search parameters were adjusted to limit the search to the common Hymenoptera-Aculeata sub-clade
(as much as possible). The resulting multi-sequence file contains over 2500 sequences.

NOTE:  Alternative "bar code" conserved sequences were also searched for via the NCBI CDD database, but
       the available data could not be downloaded or processed as easily and this approach was left for
       another day.
       
Storm Plots:
 
   To look for conserved regions in the mtDNA sequences, storm plots were created for each pair of species
-specific FASTA files.  The source python code can be found in "storm_plot.py".  The program usage is:

          python storm_plot.py Pol_d_mtDNA_CO1_cds.fasta S_invicta_mtDNA_CO1_cds.fasta
          
          python storm_plot.py Trypoxylon_politum_CO1_cds.fasta S_invicta_mtDNA_CO1_cds.fasta
          
          python storm_plot.py Trypoxylon_politum_CO1_cds.fasta Pol_d_mtDNA_CO1_cds.fasta

 
Multiple Sequence Alignment:

   Sequence alignments were run against the protein BLAST mFasta file with a S. invicta COI protein
sequence using the Needleman-Wunseh scoring matrix algorithm.  Due to large number of sequences, the 
program generates a random sample from the sequences file for processing.  Each pairwise alignment score
is recorded and used to generate overlapping histograms for each genera.  The first-order statistics
(mean and standard deviation estimates) are also generated.  

   The processing source code can be found in "align_fast.py".  The program usage is:
   
      python align_fasta.py Aculeata_COX1_protein_2.fasta S_invicta_protein_CO1.fasta


RESULTS:

All plots can be found in the 'plots' subdirectory as png images.

The plot between the ant and eusocial wasp do indicate a partial conserved region and so does the plot 
between the two wasp sequences.  No indication of a conserved sequence appears in the ant versus
mud dauber sequences.

The first-order statistics for sequence alignment scores with a S.invicta query sequence are:

Solenopsis  mean =  877.52  stdev =  268.77
Polistes    mean =  677.12  stdev =  380.21
Trypoxylon  mean =  527.34  stdev =  76.83

The combined histogram plots are consistent with these stats, but show a very broad scattering for the
Polistes genera.  

CONCLUSIONS:

The storm plots are interesting and focusing on the sub-sequence regions that yield the linear patterns
are a good place to start further investigation.

While alignment score statistics are consistent with the hypothesis, Trypoxylon species are under-
representated in the data and the means are all within two standard deviations of each other.

The results do show a stronger relationship between a eusocal wasp genera and Solenopis ants, 
but are not conclusive.  


The ground has already covered in detail here:

http://www.cell.com/current-biology/fulltext/S0960-9822(17)30325-1











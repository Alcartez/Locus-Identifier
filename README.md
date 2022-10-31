# Locus Finder

## Why ?

Developed as a part of ELAN pipeline. For finding pre-insertion locii , of around 200 nucleotide length. 

## Requirements

* Blast + 
* Biopython
* Any OS

## Installation 

1. Download the github repo and run the Locus_Identifier_script.py

## Output 

It will generate 2 output files

* output.csv will contain the query matches using nblast
* LocusIdentifier_output.csv will contain the locuses.

## Algorithm and Interpretation 

In order to interpret the results , one must understand how it works.

1. First the algorithm runs blastn using a local database developed by the algorithm itself , to find transposons (in this case Alu elements).
2. The algorithm then finds the location of the specified transposon in the input sequence(s).
3. The algorithm then finds the starting and the ending positions of the transposons. 
4. Then the algorithm proceeds to find the locus. The locus is 100 nucleotides downstream and 100 nucleotides upstream of the transposon. The locus in simple terms is the 200 nucleotide length sequence that came before the transposon insertion. 
5. Once the algorithm finds the locii , it creates a dataframe with the downstream position , upstream position , name of the transposon that matched and the locus itself. 

# PlasmidBot
A pipeline to identify plasmids in short read sequences

Plasmidbot is a pipeline I created to give a greater level of confidence of the existence of plasmid sequences when using short read sequences for DNA analysis of bacterial whole genomes.

Plasmidbot is installed using Conda to create the necessary environment and build the required dependencies.
A requirement of the pipeline is the PLSDB database among a couple of others.

The pipeline makes use of various databases and tools to search for hits to known plasmid sequences.
It then downloads the full sequences for the hits from the NCBI data base and runs an alignment.
From that alignment it determines the 'percentage coverage' against the plasmid sequence, providing a text report of the findings.

The script still requires a bit of work, especially around the 'concensus' FNA which shows matches of sequences and their particular gene names. 

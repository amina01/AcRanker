Desktop version of AcRanker

AcRanker is a machine learning system developed in python that ranks proteins in a proteome as per their Anti-CRISPR tendencies predicted using sequence features.

The method takes as input a proteome in FASTA format and returns a ranked list as per the expected Acr behavior in csv format.

Use the following command to generate predictions:
python acranker.py <input file name> <output file name>

For example,
python acranker.py input_proteome.fasta results

Following packages are required:
Biopython
sklearn
scipy


Desktop version of AcRanker

Huang, Le, Bowen Yang, Haidong Yi, Amina Asif, Jiawei Wang, Trevor Lithgow, Han Zhang, Fayyaz ul Amir Afsar Minhas, and Yanbin Yin. “AcrDB: A Database of Anti-CRISPR Operons in Prokaryotes and Viruses.” Nucleic Acids Research. Accessed October 24, 2020. https://doi.org/10.1093/nar/gkaa857.

Eitzinger, Simon, Amina Asif, Kyle E. Watters, Anthony T. Iavarone, Gavin J. Knott, Jennifer A. Doudna, and Fayyaz ul Amir Afsar Minhas. “Machine Learning Predicts New Anti-CRISPR Proteins.” Nucleic Acids Research 48, no. 9 (May 21, 2020): 4698–4708. https://doi.org/10.1093/nar/gkaa219.

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


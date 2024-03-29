עמרי נהור 211359054
ניצן רפאל מגדל 206763351
יעל בורוכוב 323029637


Bioinformatics Final Project
This project analyzes a GenBank file for Bacillus clausii and compares it to UniProt data.

Files
main.py - Main script to analyze GenBank file
part_b/part_b.py - Script to analyze UniProt data
part_c/part_c.py - Script to retrieve sequences from NCBI using accession numbers
Usage
The main scripts can be run as:

part a:
python main.py
part b:
python part_b/part_b.py 
part c:
python part_c/one.py
python part_c/two.py
python part_c/dnds.py



main.py will output statistics and plots based on analyzing the GenBank file.

part_b.py will analyze the UniProt data and output hydrophobicity statistics.

part_c is divided into 3 parts:
one.py - prints according to the standard table all the codons and their options to be synonymous and non-synonymous
two.py-checks whether there are shared genes and prints them (there are no different genes between the 2 files)
dnds.py-performs dnds calculation and prints 5 common genes.
2 attachments: April.gb February.gb that we downloaded about the corona virus.

Results
Key results from the analysis include:

Statistics on gene features in GenBank file
Comparison of genes on positive vs negative strands
Histogram of gene lengths
Analysis of TC count for each gene
Comparison of GenBank vs UniProt gene names
Hydrophobicity analysis of UniProt proteins
Contributing
Contributions to expand the analysis are welcome!

License
This project is unlicensed.

# Multiparent percentages

usage: mpp_percentages.py [-h] --variants VARIANTS --progeny-substring PROGENY_SUBSTRING --parent-substrings PARENT_SUBSTRINGS 
                               [--subset-size SUBSET_SIZE] [--output OUTPUT] [--multiallelic] 

This script was designed to find the percent representation of parents in a set of recombinant progeny based on unique variants

optional arguments:
-h, --help              show this help message and exit
--subset-size SUBSET_SIZE
                        Number of SNPs for subset bar plot.
--output OUTPUT         File name for bar plot.
--multiallelic          Indicator for biallelic variants or not.

required arguments:
    --variants VARIANTS The name of the input file (CSV). File structure: R qtl2 genotype format
    --progeny-substring PROGENY_SUBSTRING
                The substring present in all progeny
    --parent-substrings PARENT_SUBSTRINGS
                The substrings present in each parent -- comma separated 

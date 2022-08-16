import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt
from collections import OrderedDict


def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description="This script was designed to find the percent representation of " +
                                                 "parents in a set of recombinant progeny based on unique variants\n\n")

    required_named = parser.add_argument_group('required arguments')

    required_named.add_argument(
        "--variants",
        type=str,
        required=True,
        help="The name of the input file (CSV). " +
             "File structure: R qtl2 genotype format",
        action="store")

    required_named.add_argument(
        "--progeny-substring",
        type=str,
        required=True,
        help="The substring present in all progeny",
        action="store")

    required_named.add_argument(
        "--parent-substrings",
        type=str,
        required=True,
        help="The substrings present in each parent -- comma separated",
        action="store")

    parser.add_argument(
        "--subset-size",
        type=int,
        required=False,
        default=0,
        help="Number of SNPs for subset bar plot.",
        action="store")

    parser.add_argument(
        "--output",
        type=str,
        required=False,
        default="barplot.png",
        help="File name for bar plot.",
        action="store")

    parser.add_argument(
        "--multiallelic",
        required=False,
        help="Indicator for biallelic variants or not.",
        action="store_true")

    return parser.parse_args()


def get_samples_with_substring(df, substring):
    return df.loc[df.index.str.contains(substring), :]


def get_variant_set_consensus(df):
    consensus = pd.to_numeric(df.apply(pd.Series.mode).loc[0, ], downcast='integer')
    return consensus


def uniq_parent_id(series):
    """Determine most common genotype (A/B) among parents by mode, return uniq parent"""
    mode = pd.Series.mode(series).sample() # sample in case of equal bialleles
    uniq_variants.append(int(mode))
    return list(series.where(series != int(mode)).dropna().index)[0]


def sample_x_variants(df, substring, number):
    """Return sample of variants of size x for parent groups"""
    parent_df = df.loc[df.Parent == substring, :]
    if parent_df.shape[1] < number:
        number = parent_df.shape[1]
        return df.loc[df.Parent == substring, :].sample(n=number)
    else:
        return df.loc[df.Parent == substring, :].sample(n=number)


def write_subset(snp_subset):
    with open("even_snp_set.txt", "w") as f:
        for snp in snp_subset:
            f.write(snp + "\n")


def get_parent_percentages(uniq_snp_to_parent_df, parent_labels):
    parent_percentages = pd.DataFrame(columns=parent_labels)
    fast_snp_parent_dictionary = OrderedDict()
    for variant_name in list(uniq_snp_to_parent_df.index):
        series = uniq_snp_to_parent_df.loc[variant_name, ]
        fast_snp_parent_dictionary[variant_name] = (series.Parent, series.Variant)
    for mb in tqdm(mb_progeny_uniq.index):
        mb_percentage = pd.Series([0] * len(parent_labels), name=mb, index=parent_labels)
        for variant_name in list(uniq_snp_to_parent_df.index):
            uniq_variant = fast_snp_parent_dictionary[variant_name]
            if uniq_variant[1] == mb_progeny_uniq.loc[mb, variant_name]:
                mb_percentage[parent_labels.index(uniq_variant[0])] += 1
        parent_percentages = parent_percentages.append(mb_percentage/sum(mb_percentage))
    return parent_percentages


def plot_percentages_sorted_by_parent(df, sort_parent, output):
    """Given dataframe of percentage and parent to sort by, generate figure"""
    df.sort_values(by=sort_parent).plot(kind="bar", stacked=True, edgecolor="none", width=1.0)
    plt.tick_params(
        axis='x',           # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        bottom=False,       # ticks along the bottom edge are off
        top=False,          # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    plt.savefig(output)


if __name__ == '__main__':
    arguments = parse_arguments()
    if not any([arguments.output.endswith(i) for i in ["eps", "pdf", "pgf", "png", "ps", "raw", "rgba", "svg", "svgz"]]):
        sys.stderr.write("Output must end with EPS|PDF|PNG|etc.\n")
        sys.exit(1)
    sys.stderr.write("Reading file\n")
    variants = pd.read_csv(arguments.variants, index_col=0, comment="%")
    print(variants.head())
    if arguments.multiallelic:
        recode_allele = {"A": 0, "C": 1, "G": 2, "T": 3, 
                "-": -1, "NA": -1, np.nan: -1, "N": -1, 
                "K":4, "Y": 5, "R":6, "B":7, "W":8, "M": 9, "S": 10}
    else:
        recode_allele = {"A": 0, "B": 1, "-": -1, "NA": -1, np.nan: -1}
    sys.stderr.write("Recoding variants\n")
    variants = variants.applymap(lambda x: recode_allele[x])
    progeny_substring = arguments.progeny_substring
    parents = arguments.parent_substrings.split(",")

    sys.stderr.write("Splitting sample sets\n")
    sample_sets = []
    parent_consensus_sets = []
    for parent in parents:
        sample_set = get_samples_with_substring(variants, parent)
        sample_sets.append(sample_set)
        parent_consensus = get_variant_set_consensus(sample_set)
        parent_consensus_sets.append(parent_consensus)
    sample_sets.append(get_samples_with_substring(variants, progeny_substring))

    sys.stderr.write("Isolating parent specific SNPs\n")
    parent_genotypes = pd.DataFrame(parent_consensus_sets, index=parents)
    parent_genotypes_sums = parent_genotypes.sum(axis=0)
    uniq_parent_snps = parent_genotypes_sums[(parent_genotypes_sums == 1) | (parent_genotypes_sums == len(parents) - 1)]
    uniq_parent_genotypes = parent_genotypes.loc[:, list(uniq_parent_snps.index)]

    uniq_variants = []

    uniq_snp_to_parent = uniq_parent_genotypes.apply(uniq_parent_id)
    uniq_snp_to_parent.value_counts()

    uniq_snp_to_parent_encoding = pd.DataFrame(uniq_snp_to_parent, columns=["Parent"])
    uniq_snp_to_parent_encoding["Variant"] = uniq_variants

    mb_progeny_uniq = sample_sets[-1].loc[:, list(uniq_snp_to_parent_encoding.index)]

    percentages = get_parent_percentages(uniq_snp_to_parent_encoding, parents)
    plot_percentages_sorted_by_parent(percentages, parents[0], arguments.output)

    if arguments.subset_size >= 1:

        parent_sample_subsets = []
        for parent in parents:
            parent_sample_subsets.append(sample_x_variants(uniq_snp_to_parent_encoding, parent, arguments.subset_size))
    
        even_snp_set = []
        for i in parent_sample_subsets:
            even_snp_set = even_snp_set + list(i.index)
        write_subset(even_snp_set)
    
        
        even_percentages = get_parent_percentages(uniq_snp_to_parent_encoding.loc[even_snp_set, ], parents)
        plot_percentages_sorted_by_parent(even_percentages, parents[0], arguments.output.replace(".png", ".even.png"))
        

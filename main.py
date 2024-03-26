from operator import ne, pos
from Bio import SeqIO
import matplotlib.pyplot as plt
from numpy import average, mean, negative
from typing import List
import pandas as pd
import numpy as np


def create_type_table():
    counts = {}

    # bacilus_file = open("./Bacillus clausii.gb")
    # print(bacilus_file.read
    # ())
    sequence_count = 0
    # records = SeqIO.parse(bacilus_file, "genbank")
    for record in SeqIO.parse(open("./Bacillus clausii.gb"), "genbank"):
        for feature in record.features:
            # print(feature.type)
            if feature.type in counts:
                counts[feature.type] += 1
            else:
                counts[feature.type] = 1
        # print(sequence)
        print(sequence_count)
    return counts


def create_histogram():
    gene_lengths = []
    for record in SeqIO.parse(open("./Bacillus clausii.gb"), "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                gene_lengths.append(feature.location.end - feature.location.start)
    print_histogram("gene_lengths", gene_lengths, ["gene length", "count"])


def print_histogram(name: str, values, column_names: List[str]):
    plt.hist(
        values,
        bins=range(int(min(values)), int(max(values)) + 1),
        edgecolor="black",
    )
    # Add labels and title
    plt.xlabel(column_names[0])
    plt.ylabel(column_names[1])
    plt.title(name)
    # Show plot
    plt.show()


def compare_strands():
    positive_strand = 0
    negative_strand = 0
    cds_positive_lengths = []
    cds_negative_lengths = []
    non_cds_positive = []
    non_cds_negative = []
    positive_strand_cds = 0
    negative_strand_cds = 0
    for record in SeqIO.parse(open("./Bacillus clausii.gb"), "genbank"):
        # print(record.seq)
        for feature in record.features:
            if feature.location.strand == -1:
                if feature.type == "CDS":
                    negative_strand += 1
                    negative_strand_cds += 1
                    cds_negative_lengths.append(
                        len(
                            record.seq[
                                feature.location.start : feature.location.end
                            ].reverse_complement()
                        )
                    )
                elif feature.type == "gene":
                    negative_strand += 1
                    non_cds_negative.append(
                        len(
                            record.seq[
                                feature.location.start : feature.location.end
                            ].reverse_complement()
                        )
                    )

            elif feature.location.strand == 1:
                if feature.type == "CDS":
                    positive_strand += 1
                    positive_strand_cds += 1
                    cds_positive_lengths.append(
                        feature.location.end - feature.location.start
                    )
                elif feature.type == "gene":
                    positive_strand += 1
                    non_cds_positive.append(
                        feature.location.end - feature.location.start
                    )

    print(
        "positive_strand: ",
        positive_strand,
        "negative_strand: ",
        negative_strand,
        "positive_strand_with_cds: ",
        positive_strand_cds,
        "negative_strand_cds: ",
        negative_strand_cds,
    )
    # histograms
    print_histogram(
        "positive_cds", cds_positive_lengths, ["length", "precentege"]
    )  # cds positive histogram.
    print_histogram(
        "negative_cds", cds_negative_lengths, ["length", "precentege"]
    )  # cds negative histogram.
    # print(non_cds_positive)
    print_histogram(
        "positive_non_cds", non_cds_positive, ["length", "precentege"]
    )  # non cds positive
    print_histogram(
        "negative_non_cds", non_cds_negative, ["length", "precentege"]
    )  # non cds negative
    # this is to print the stats for the cds statitatics.
    print_statistic_table(cds_positive_lengths, cds_negative_lengths)
    # this prints the statistics of non cds genes.
    print_statistic_table(non_cds_positive, non_cds_negative)


def print_statistic_table(positive, negative):
    # min max avarage and standart deviasion
    min = [np.min(positive), np.min(negative)]
    max = [np.max(positive), np.max(negative)]
    average = [np.mean(positive), np.mean(negative)]
    std = [np.std(positive), np.std(negative)]
    # Create a DataFrame

    # Calculate statistics
    statistics = {
        "min": min,
        "max": max,
        "average": average,
        "standard deviation": std,
    }

    # Create a new DataFrame with statistics.
    statistics_df = pd.DataFrame(statistics)
    statistics_df.rows = ["positive", "negative"]
    print(statistics_df)
    #


# this function gets the percent of tc per percentage of each gene
def get_whole_dna_strand():
    tc_count_in_record = 0
    record_len = 0
    tc_proteins = pd.DataFrame(columns="gene,strand,start,end,tc_count".split(","))
    for record in SeqIO.parse(open("./Bacillus clausii.gb"), "genbank"):
        seq_record = record.seq
        tc_count_in_record = count_tc(seq_record)
        record_len = len(record)
        for feature in record.features:
            if feature.type == "CDS":
                # print(feature.qualifiers)
                if "gene" in feature.qualifiers:
                    count = 0
                    if feature.location.strand == 1:
                        count = (
                            count_tc(
                                record.seq[
                                    feature.location.start : feature.location.end
                                ]
                            )
                            * 100
                            / (feature.location.end - feature.location.start)
                        )

                    elif feature.location.strand == -1:
                        count = (
                            count_tc(
                                record.seq[
                                    feature.location.start : feature.location.end
                                ].reverse_complement()
                            )
                            * 100
                            / (feature.location.end - feature.location.start)
                        )

                    data_to_append = pd.Series(
                        {
                            "gene": feature.qualifiers["gene"][0],
                            "strand": feature.location.strand,
                            "start": feature.location.start,
                            "end": feature.location.end,
                            "tc_count": count,
                        }
                    )
                    tc_proteins = pd.concat([tc_proteins, data_to_append.to_frame().T])
        print(tc_proteins)
    print("tc_count_in_record", tc_count_in_record / record_len * 100)
    print("tc_protein", np.mean(tc_proteins["tc_count"]))
    print("tc_protein", tc_proteins)
    print_histogram(
        "tc_proteins", tc_proteins["tc_count"], ["precentege", "count per percentage"]
    )  # cds positive histogram."tc_proteins)
    # Sort genes by average values while keeping track of original indices
    sorted_tc_proteins = pd.DataFrame(tc_proteins).sort_values(by="tc_count")
    print(sorted_tc_proteins[-6:-1])
    print(sorted_tc_proteins[:5])
    return tc_proteins


def count_tc(seq_record: str):
    tc_count_in_record = seq_record.count("T") + seq_record.count("C")
    return tc_count_in_record


print(create_type_table())
create_histogram()
compare_strands()
get_whole_dna_strand()

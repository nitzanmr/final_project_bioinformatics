from Bio import SeqIO
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from typing import List


hydrophobic_aa = ["A", "V", "L", "I", "M", "F", "W", "P"]


def count_tc(seq_record: str):
    tc_count_in_record = seq_record.count("T") + seq_record.count("C")
    return tc_count_in_record


def get_whole_dna_strand():
    tc_count_in_record = 0
    record_len = 0
    tc_proteins = pd.DataFrame(columns="gene,strand,start,end,tc_count".split(","))
    for record in SeqIO.parse(open("../Bacillus clausii.gb"), "genbank"):
        seq_record = record.seq
        tc_count_in_record = count_tc(seq_record)
        record_len = len(record)
        for feature in record.features:
            if feature.type == "CDS":
                # print(feature.qualifiers)
                if "gene" in feature.qualifiers:
                    count = (
                        count_tc(record[feature.location.start : feature.location.end])
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
    # print_histogram(
    #     "tc_proteins", tc_proteins["tc_count"], ["precentege", "count per percentage"]
    # )  # cds positive histogram."tc_proteins)
    # Sort genes by average values while keeping track of original indices
    sorted_tc_proteins = pd.DataFrame(tc_proteins).sort_values(by="tc_count")
    print(sorted_tc_proteins[-6:-1])
    print(sorted_tc_proteins[:5])
    return tc_proteins


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


def print_statistic_table_protiens(proteins_lengths: List[int]):
    # min max avarage and standart deviasion
    # Create a DataFrame
    min = np.min(proteins_lengths)
    max = np.max(proteins_lengths)
    average = np.mean(proteins_lengths)
    std = np.std(proteins_lengths)
    # Calculate statistics
    statistics = {
        "min": min,
        "max": max,
        "average": average,
        "standard deviation": std,
    }

    # Create a new DataFrame with statistics.
    print(statistics)
    #


# this is an analysis for the uniprot file of bacilluis.
def uniprot_bacillus(path: str):
    records = []
    records_len = []
    for record in SeqIO.parse(path, "fasta"):
        print(record.seq)
        records_len.append(len(record.seq))
        pattern = r"GN=([^ ]+)"
        match = re.search(pattern, record.description)
        if match:
            gene_name = match.group(1)  # Extract the gene name from the matched group
            # print("Gene Name:", gene_name)
            records.append(gene_name)
            # else:
            # print("Gene Name not found.")
    print_histogram("records_len", records_len, ["precentege", "count per percentage"])
    # print(records)
    return records


def check_hedrophobic_genes(path: str):
    hydrophobic_count_arr = []
    i = 0
    for record in SeqIO.parse(open(path), "fasta"):
        count = 0
        for aa in record:
            if aa in hydrophobic_aa:
                count += 1
        hydrophobic_count_arr.append(count / len(record) * 100)
        i += 1
    print(hydrophobic_count_arr)
    print_statistic_table_protiens(hydrophobic_count_arr)
    return hydrophobic_count_arr


def compare_ids(uniprot_ids, genbank_ids):
    # genbank_ids["gene"] = genbank_ids["gene"].str[0].str.upper()
    # print(genbank_ids["gene"].isin(uniprot_ids).values)
    # print(genbank_ids)
    genbank_ids["gene"] = genbank_ids["gene"].mask(
        genbank_ids["gene"].isin(uniprot_ids)
    )
    genbank_ids = genbank_ids.dropna(subset=["gene"])

    print(genbank_ids)
    return genbank_ids


ids = get_whole_dna_strand()
print(ids)
uniprot_ids = uniprot_bacillus("./uniprotkb_bacillus_clausii_2024_03_13.fasta")
genes_not_in_uniprot = compare_ids(
    uniprot_ids, ids.copy()
)  # there are less genes in genbank not in uniprot

check_hedrophobic_genes("./uniprotkb_bacillus_clausii_2024_03_13.fasta")
# print(ids["gene"])
# uniprot_not_in_genbank = [x for x in uniprot_ids if x not in ids["gene"].values]
# print(len(uniprot_ids))
# print(uniprot_not_in_genbank)  # there are a lot more proteins in uniprot not in genbank

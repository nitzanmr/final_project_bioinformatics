from Bio import Entrez, SeqIO, Align
from typing import List
from icecream import ic
from matplotlib.pyplot import axis
import pandas as pd
import re
import random
from Bio.codonalign import CodonAlignment
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from Bio.Data.CodonTable import generic_by_id
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

gencode = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "_",
    "TAG": "_",
    "TGC": "C",
    "TGT": "C",
    "TGA": "_",
    "TGG": "W",
}


def get_fasta_files(accesion_numbers: List[str]):
    Entrez.email = "nitzanmr@gmail.com"
    fasta_files = []
    returned_sequences = []
    for i in accesion_numbers:
        fasta_files.append(
            Entrez.efetch(db="nucleotide", id=i, rettype="genbank", retmode="text")
        )
        # for record in SeqIO.parse(fasta_files[-1], "fasta"):
        #     # print("%s %i" % (record.seq, len(record)))
        #     returned_sequences.append(record)
    return fasta_files


def translate(seq):
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            protein += gencode[codon]
    return protein


#
def dnds_to_genes(gene1: str, gene2: str, strand="+"):
    # Perform nucleotide sequence alignment
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(str(gene1), str(gene2), strand="+")
    aligned_seq1, aligned_seq2 = alignments[0]
    a = SeqRecord(CodonSeq(str(aligned_seq1)), id="Alpha")
    b = SeqRecord(CodonSeq(str(aligned_seq2)), id="Beta")
    alignment = CodonAlignment([a, b])


# print(alignment)
# # aligned_seq1 = translate(gene1)
# # aligned_seq2 = translate(gene2)
# print(aligned_seq1.__len__(), aligned_seq2.__len__())
# # translated_seq = Seq(seq1).translate(table=gencode, to_stop=True)
# # Extract the aligned sequences
# dn, ds = cal_dn_ds(aligned_seq1, aligned_seq2, codon_table=None)
# return dn, ds
# aligner = Align.PairwiseAligner()
# aligner.align(gene1, gene2)
# seq1 = CodonSeq(gene1)
# translated_seq = Seq(seq1).translate(table=gencode, to_stop=True)
# print(translated_seq)
# print(seq1)
# seq2 = CodonSeq(gene2)
# print(generic_by_id[16])
# dN, dS = cal_dn_ds(seq1, seq2, codon_table=generic_by_id[11])
# dN_dS_ratio = float(dN / dS)
# print("dN:%0.3f " % dN)
# print("dS:%0.3f " % dS)
# print("dN/dS:%0.3f " % dN_dS_ratio)


fasta = ["MZ054892.1", "PP348372.1"]
fasta_files = get_fasta_files(fasta)
# icecream.ic(fasta_files)
# print(fasta_files)
stats_data = []
columns = ["length", "number_of_genes", "number_of_proteins", "gene_names", "gene_seq"]
stats = pd.DataFrame(columns=columns)
for fasta_file in fasta_files:
    number_of_proteins = 0
    number_of_genes = 0
    gene_names = []
    protein_names = []
    gene_seq = []
    for record in SeqIO.parse(fasta_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" or feature.type == "gene":
                if "gene" in feature.qualifiers:
                    if (
                        feature.qualifiers["gene"][0] not in gene_names
                        and feature.qualifiers["gene"][0] != "S"
                        and feature.qualifiers["gene"][0] != "ORF1ab"
                    ):
                        gene_names.append(feature.qualifiers["gene"][0])

                        # Extract gene sequence based on strand
                        if feature.location.strand == 1:
                            gene_seq.append(
                                record.seq[
                                    feature.location.start : feature.location.end
                                ]
                            )
                        elif feature.location.strand == -1:
                            gene_seq.append(
                                record.seq[
                                    feature.location.start : feature.location.end
                                ].reverse_complement()
                            )

                        if feature.type == "CDS":
                            number_of_proteins += 1
                        else:
                            number_of_genes += 1
        stats_data.append(
            [len(record), number_of_genes, number_of_proteins, gene_names, gene_seq]
        )

        stats = pd.DataFrame(stats_data, columns=columns)
        # print(gene_names, protein_names)
        # print(number_of_proteins, number_of_genes)

        # print(stats)
simaller_genes = 0
simaller_genes_names = []
for i in stats["gene_names"][0]:
    if i in stats["gene_names"][1]:
        simaller_genes += 1
        simaller_genes_names.append(i)
picked_genes = random.sample(simaller_genes_names, 5)
pd_picked_genes = pd.DataFrame(picked_genes, columns=["gene_names"])
gene_names = stats["gene_names"]
gene_seq = stats["gene_seq"]

new_df = []
for i in range(len(gene_names)):
    combined = []
    for j in range(len(gene_names[i])):
        # print(gene_names[i][j], gene_seq[i][j])
        combined.append([gene_names[i][j], gene_seq[i][j]])
    new_df.append(pd.DataFrame(combined, columns=["gene_names", "gene_seq"]))
    new_df[-1] = new_df[-1].merge(
        right=pd_picked_genes, right_on="gene_names", left_on="gene_names"
    )
print(new_df)
for i in range(5):
    print(len(new_df[0]["gene_seq"][i]), len(new_df[1]["gene_seq"][i]))
    dnds_to_genes(new_df[0]["gene_seq"][i], new_df[1]["gene_seq"][i])
# records_len.append(len(record.seq))
# pattern = r"N=([^ ]+)"
# match = re.search(pattern, record.description)
# if match:
#     gene_name = match.group(1)  # Extract the gene name from the matched group
#     print(gene_name)
#         stats["number_of_genes"] += 1
#     if record.type == "CDS":
#         stats["number_of_proteins"] += 1
# for record in fasta_files[1]:
#     stats["length"] = len(record)
#     if record.type == "gene":
#         stats["number_of_genes"] += 1
#     if record.type == "CDS":
#         stats["number_of_proteins"] += 1

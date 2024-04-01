# Importing necessary modules and functions
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio import Align
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from Bio.Data import CodonTable
from Bio.Align import PairwiseAligner
import pandas as pd


# Function to load sequences from GenBank files
def load_genbank_sequences(april_file, february_file):
    april_record = SeqIO.read(april_file, "genbank")
    february_record = SeqIO.read(february_file, "genbank")
    return april_record, february_record

# Function to retrieve sequence of a specific gene
def get_gene_sequence(record, gene_name):
    for feature in record.features:
        if feature.type == "CDS" and feature.qualifiers.get("gene", [""])[0] == gene_name:
            if feature.location.strand == -1:
                return record.seq[feature.location.start:feature.location.end].reverse_complement()
            return feature.location.extract(record.seq)
    return None

# Function to translate a nucleotide sequence into a protein sequence
def translate_sequence(seq):
    return str(seq.translate(table="Standard"))

# Function to align two protein sequences
def align_protein_sequences(prot_seq1, prot_seq2):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignments = aligner.align(prot_seq1, prot_seq2)
    return alignments[0]

# Function to convert protein sequence back to codon sequence based on DNA sequence
def protein_to_codons(protein_sequence, dna_sequence):
    dna_seq = ""
    dna_seq_index = 0
    stop_codons = {"TAA", "TAG", "TGA"}  # Stop codons
    
    for amino_acid in protein_sequence:
        if amino_acid == '-':
            dna_seq += '---'  # If "-" is encountered, add "---" to the DNA sequence
        else:
            codon = dna_sequence[dna_seq_index:dna_seq_index + 3]
            if codon in stop_codons:
                break  # Stop the loop if a stop codon is encountered
            dna_seq += codon
            dna_seq_index += 3  # Advance the position in dna_sequence by 3 characters
            
    return dna_seq

# Function to calculate dN and dS values between two codon sequences
def calculate_dn_ds(seq1, seq2):
    seq1 = CodonSeq(seq1)
    seq2 = CodonSeq(seq2)
    dn, ds = cal_dn_ds(seq1, seq2)
    return dn, ds

# Function to determine the type of selection based on dN and dS values
def determine_selection_type(dn, ds):
    if dn == 0 and ds == 0:
        return "Neutral"
    elif ds == 0:
        return "Positive"
    else:
        return "Negative" if dn / ds < 1 else "Positive"

# Function to retrieve details of a gene
def get_gene_details(record, gene_name):
    gene_details = {}
    for feature in record.features:
        if feature.type == "CDS" and feature.qualifiers.get("gene", [""])[0] == gene_name:
            gene_details["Role"] = feature.qualifiers.get("product", ["Unknown"])[0]
            gene_details["Start"] = feature.location.start
            gene_details["End"] = feature.location.end
            return gene_details
    return {"Role": "Unknown", "Start": "-", "End": "-"}

# Function to analyze genes and return a DataFrame with gene data
def analyze_genes(list_of_genes, april_file, february_file):
    gene_data = []
    april_record, february_record = load_genbank_sequences(april_file, february_file)
    
    for gene_name in list_of_genes:
        seq1 = get_gene_sequence(april_record, gene_name)
        seq2 = get_gene_sequence(february_record, gene_name)
        if seq1 is None or seq2 is None:
            continue  # Skip if sequence not found
        
        length_april = len(seq1)
        length_february = len(seq2)

        translated_seq1 = translate_sequence(seq1)
        translated_seq2 = translate_sequence(seq2)
        
        alignments = align_protein_sequences(translated_seq1, translated_seq2)
        alignment = alignments  # considering the first alignment
        codon_sequence1 = protein_to_codons(alignment[0], seq1)
        codon_sequence2 = protein_to_codons(alignment[1], seq2)

        dn, ds = calculate_dn_ds(codon_sequence1, codon_sequence2)
        dn_ds_selection = determine_selection_type(dn, ds)
        
        gene_details_april = get_gene_details(april_record, gene_name)
        gene_details_february = get_gene_details(february_record, gene_name)
        
        gene_data.append({
            "Gene Name": gene_name,
            "Role": gene_details_april["Role"],
            "Start (April)": gene_details_april["Start"],
            "End (April)": gene_details_april["End"],
            "Start (February)": gene_details_february["Start"],
            "End (February)": gene_details_february["End"],
            "Length in April": length_april,
            "Length in February": length_february,
            "dN/dS Selection": dn_ds_selection
            
        })

    gene_df = pd.DataFrame(gene_data)
    return gene_df

# Main function to execute the analysis
def main():
    april_file = "April.gb"
    february_file = "February.gb"
    # list_of_genes = ["ORF6","ORF3a","ORF10","ORF7b","E","ORF8","N","M","ORF7a"]
    list_of_genes = ["M","N","ORF3a","ORF7a","ORF10"]
    
    gene_table = analyze_genes(list_of_genes, april_file, february_file)
    print(gene_table)

if __name__ == "__main__":
    main()

from Bio import SeqIO

# Function to get genome length
def get_genome_length(file):
    for record in SeqIO.parse(file, "genbank"):
        return len(record.seq)

# Function to get number of genes and protein coding genes
def get_gene_info(file):
    total_genes = 0
    protein_coding_genes = 0
    gene_names = set()

    # print("hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh")
    for record in SeqIO.parse(file, "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                total_genes += 1
                if 'gene' in feature.qualifiers:
                    gene_names.add(feature.qualifiers['gene'][0])
                    # print(feature.qualifiers['gene'][0],len(feature.location))
            elif feature.type == "CDS":
                protein_coding_genes += 1
                if 'gene' in feature.qualifiers:
                    gene_names.add(feature.qualifiers['gene'][0])


    return total_genes, protein_coding_genes, gene_names

# Function to compare genomes and find common and unique genes
def compare_genomes(april_file, february_file):
    april_genes = set()
    february_genes = set()

    # Get gene names for April genome
    for record in SeqIO.parse(april_file, "genbank"):
        for feature in record.features:
            if feature.type in ["gene", "CDS"]:
                if 'gene' in feature.qualifiers:
                    april_genes.add(feature.qualifiers['gene'][0])

    # Get gene names for February genome
    for record in SeqIO.parse(february_file, "genbank"):
        for feature in record.features:
            if feature.type in ["gene", "CDS"]:
                if 'gene' in feature.qualifiers:
                    february_genes.add(feature.qualifiers['gene'][0])

    common_genes = april_genes.intersection(february_genes)
    unique_april_genes = april_genes - february_genes
    unique_february_genes = february_genes - april_genes

    return len(common_genes), list(unique_april_genes), list(unique_february_genes)

# Main function
def main(april_file, february_file):
    # A. Get genome length
    april_genome_length = get_genome_length(april_file)
    february_genome_length = get_genome_length(february_file)
    print("A. Genome Length:")
    print("April Genome Length:", april_genome_length)
    print("February Genome Length:", february_genome_length)

    # B. Get number of genes and protein coding genes
    april_total_genes, april_protein_coding_genes, april_gene_names = get_gene_info(april_file)
    february_total_genes, february_protein_coding_genes, february_gene_names = get_gene_info(february_file)
    print("\nB. Gene Information:")
    print("April Total Genes:", april_total_genes)
    print("April Protein Coding Genes:", april_protein_coding_genes)
    print("April Gene Names:", ", ".join(april_gene_names))
    print("February Total Genes:", february_total_genes)
    print("February Protein Coding Genes:", february_protein_coding_genes)
    print("February Gene Names:", ", ".join(february_gene_names))

    # C. Compare genomes
    common_genes, unique_april_genes, unique_february_genes = compare_genomes(april_file, february_file)
    print("\nC. Comparison of Genomes:")
    print("Number of Common Genes:", common_genes)
    print("Genes unique to April Genome:", unique_april_genes)
    print("Genes unique to February Genome:", unique_february_genes)

if __name__ == "__main__":
    april_file = "April.gb"
    february_file = "February.gb"
    main(april_file, february_file)

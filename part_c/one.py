from Bio.Data import CodonTable
# print(Bio.Data.CodonTable.standard_dna_table)

def count_mutations_with_same_aa(codon, codon_table):
    count = 0
    count_stop=0
    original_aa = codon_table.forward_table[codon]
    for i in range(len(codon)):
        for base in "TGCA":
            mutated_codon = codon[:i] + base + codon[i+1:]
            mutated_aa = codon_table.forward_table.get(mutated_codon, None)
            if None==mutated_aa:
                count_stop+=1
            if mutated_codon != codon and mutated_aa == original_aa:
                count += 1
    return count,count_stop

# Get the standard DNA codon table
standard_table = CodonTable.standard_dna_table

# Iterate over the codons and count mutations with the same amino acid
for codon in standard_table.forward_table.keys():
    Synonyms,count_stop = count_mutations_with_same_aa(codon, standard_table)
    # print(f"Codon: {codon}, Mutations with same amino acid: {mutation_count}, Stop: {count_stop}")
    not_synonyms=(9-count_stop-Synonyms)
    percent_synonyms=(3*Synonyms)/(9-count_stop)
    percent_not_synonyms=(3*(not_synonyms)/(9-count_stop))
    print(f"Codon: {codon}, Percentage of synonyms: {percent_synonyms}, Percentage of not synonyms: {percent_not_synonyms}")
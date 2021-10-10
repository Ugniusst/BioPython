
from GenomeComparison import GenomeComparison
from ProteinList import ProteinList

genome_comparison = GenomeComparison()
protein_lister = ProteinList()

protein_list = protein_lister.get_all_proteins_list()

file_names = ["bacterial1.fasta","bacterial2.fasta","bacterial3.fasta","bacterial4.fasta", "mamalian1.fasta","mamalian2.fasta","mamalian3.fasta","mamalian4.fasta"]

for x in range(len(file_names)):
    for y in range(x+1, len(file_names)):
        genome_comparison.print_distances_matrix(file_names[x], file_names[y], protein_list)

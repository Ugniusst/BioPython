from FastaParser import FastaParser
from FrequencyCounter import FrequencyCounter
from OpenReaderFrameFinder import ORFsFinder
from DistanceCounter import DistanceCounter
from tabulate import tabulate

class GenomeComparison:

    fasta_parser = FastaParser()
    orf_finder = ORFsFinder()
    frequency_counter = FrequencyCounter()
    distances_counter = DistanceCounter()


    def print_distances_matrix(self, genome_file_1, genome_file_2, protein_list):

        my_seq = self.fasta_parser.parse_fasta_file(genome_file_1)
        orfs = self.orf_finder.get_all_ORFS(my_seq)
        dicodon_frequencies = self.frequency_counter.get_all_dicodon_frequencies(orfs, protein_list)

        my_seq_2 = self.fasta_parser.parse_fasta_file(genome_file_2)
        orfs_2 = self.orf_finder.get_all_ORFS(my_seq_2)
        dicodon_frequencies_2 = self.frequency_counter.get_all_dicodon_frequencies(orfs_2, protein_list)

        distances = self.distances_counter.get_all_dicodon_distances(dicodon_frequencies, dicodon_frequencies_2, protein_list)

        print(genome_file_1.split(".")[0] + " - " + genome_file_2.split(".")[0])
        print("20") 
        print(tabulate(distances, tablefmt="plain"))

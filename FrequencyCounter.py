from Bio.Seq import Seq

class FrequencyCounter:

    def concatenate_dna_sequences(self, dna_sequencies):
        concatenated = Seq("")
        for dna_sequence in dna_sequencies:
            concatenated += dna_sequence
        return concatenated

    def count_codon_frequency(self, searchable_codon, dna_sequence):
        counter = 0
        for x in range(0, len(dna_sequence),3):
            codon = dna_sequence[x:x+3]
            if codon == searchable_codon:
                counter += 1
        frequency = (counter/(len(dna_sequence)/3))*1000
        return frequency

    def count_all_codon_frequencies(self, dna_sequence):
        frequencies = [[[0 for j in range(4)] for k in range(4)] for l in range(4)]
        symbols = ["T", "C", "A", "G"]
        for s1 in range(4):
            for s2 in range(4):
                for s3 in range(4):
                    codon = symbols[s1]+symbols[s2]+symbols[s3]
                    frequencies[s1][s2][s3] = self.count_codon_frequency(Seq(codon), dna_sequence)
        return frequencies

    def count_dicodon_frequency(self, searchable_dicodon, dna_sequence):
        counter = 0
        dna_sequence = dna_sequence.translate()
        for x in range(0, len(dna_sequence)-1):
            dicodon = dna_sequence[x:x+2]
            if dicodon == searchable_dicodon:
                counter += 1
        frequency = (counter/((2 * len(dna_sequence))/2 -1))*1000
        return frequency

    def count_all_dicodon_frequencies(self, dna_sequence, proteins_list):
        dicodon_frequencies = [[0 for j in range(len(proteins_list)+1)] for k in range(len(proteins_list))]
        for p1 in range(len(proteins_list)):
            dicodon_frequencies[p1][0] = proteins_list[p1]
            for p2 in range(len(proteins_list)):
                dicodon = proteins_list[p1] + proteins_list[p2]
                dicodon_frequencies[p1][p2+1] = self.count_dicodon_frequency(dicodon, dna_sequence)
        return dicodon_frequencies

    def get_all_dicodon_frequencies(self, sequencies, protein_list):
        sequence_string = self.concatenate_dna_sequences(sequencies)
        dicodon_frequencies = self.count_all_dicodon_frequencies(sequence_string, protein_list)
        return dicodon_frequencies

    def get_all_codon_frequencies(self, sequencies):
        sequence_string = self.concatenate_dna_sequences(sequencies)
        dicodon_frequencies = self.count_all_codon_frequencies(sequence_string)
        return dicodon_frequencies
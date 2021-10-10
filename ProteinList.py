from Bio.Seq import Seq

class ProteinList:

    stop_Codons = [Seq("TAA"),Seq("TAG"),Seq("TGA")]

    def get_protein(self, codon):
        codon = codon.translate()
        return codon

    def get_all_proteins_list(self):
        proteins_list = []
        symbols = ["T", "C", "A", "G"]
        for s1 in range(4):
            for s2 in range(4):
                for s3 in range(4):
                    codon = symbols[s1]+symbols[s2]+symbols[s3]
                    if codon not in self.stop_Codons:
                        protein = self.get_protein(Seq(codon))
                        if protein not in proteins_list:
                            proteins_list.append(protein)
        return proteins_list
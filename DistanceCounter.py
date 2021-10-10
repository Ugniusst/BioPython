import math
from Bio.Seq import Seq

class DistanceCounter:

    proteins_list = []
    proteins_name_list = ["Phe", "Leu", "Ser", "Tyr", "Cys", "Trp", "Pro", "His", "Qln", "Arg", "Ile", "Met", "Thr", "Asn", "Lys", "Val", "Ala", "Asp", "Glu", "Gly"]

    def count_all_dicodon_distances(self, matrix1, matrix2):
        dicodon_distances = [[0 for j in range(len(self.proteins_list)+1)] for k in range(len(self.proteins_list))]
        for p1 in range(len(self.proteins_list)):
            dicodon_distances[p1][0] = self.proteins_name_list[p1]
            for p2 in range(len(self.proteins_list)):
                if matrix1[p1][p2+1] != 0 and matrix2[p1][p2+1] != 0:
                    fraction1 = matrix1[p1][p2+1]/(matrix1[p1][p2+1]+matrix2[p1][p2+1])
                    fraction2 = matrix2[p1][p2+1]/(matrix1[p1][p2+1]+matrix2[p1][p2+1])
                    dicodon_distances[p1][p2+1] = round(abs(fraction1-fraction2),4)
                else:
                    dicodon_distances[p1][p2+1] = 0
        return dicodon_distances

    def get_all_dicodon_distances(self, matrix1, matrix2, proteins_list):
        self.proteins_list = proteins_list
        dicodon_distances = self.count_all_dicodon_distances(matrix1, matrix2)
        return dicodon_distances

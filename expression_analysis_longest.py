#!/bin/env python
# coding: utf-8

import sys
import traceback
from collections import defaultdict

from expression_analysis import ExpressionAnalysis


class ExpressionAnalysisLongest(ExpressionAnalysis):
    def start(self):
        self.get_args()
        try:
            self.check_fastq_path()
            self.check_file_path()
            self.echo_input_params()
            self.read_trimmomatic_path()
            self._log("=== Analysis start. ===")
            if self.hisat2_index_path is None:
                self.run_hisat2_build()
            self.run_trimmomatic()
            self.run_hisat2()
            self.run_samtools_sort()
            self.run_stringtie()
            self.format_gtf_to_tsv()
            self.merge_tsv()
            self.find_longest_isoform()
            self.get_entrez_info()
            self.get_go_info()
            self.dump_tsv()
            self.dump_tsv_FPKM()
            self.dump_tsv_TPM()
            self.dump_tsv_Coverage()
            self._log("=== Analysis finish. ===")
        except:
            traceback.print_exc()
            sys.exit(1)

        return True

    def find_longest_isoform(self):
        self._log("Find longest isoform.")
        gene_ind = self.columns.index("Gene_ID")
        width_ind = self.columns.index("Width")
        chromosome_ind = self.columns.index("Chromosome")
        start_ind = self.columns.index("Start")
        d_gene_id = defaultdict(list)
        for i in range(len(self.data)):
            row = self.data[i]
            gene_id = row[gene_ind]
            width = row[width_ind]
            d_gene_id[gene_id].append([width, i])

        new_data = []
        for value in d_gene_id.values():
            value_sorted = sorted(value, key=lambda x: -x[0])
            longest_ind = value_sorted[0][1]
            new_data.append(self.data[longest_ind])

        new_data_sorted = sorted(new_data, key=lambda x: (x[chromosome_ind],
                                                          x[start_ind]))
        self.data = new_data_sorted
        self.all_refseq_id = [row[0] for row in self.data]

        return True


if __name__ == "__main__":
    my_expression_analysis = ExpressionAnalysisLongest()
    my_expression_analysis.start()


#!/bin/env python
# coding: utf-8

import sys
import traceback
import argparse
from collections import defaultdict
import mygene
import os

from expression_analysis import ExpressionAnalysis


class ExpressionAnalysisBam(ExpressionAnalysis):
    def __init__(self):
        self.l_bam_path = []
        self.gtf_path = None
        self.output_dir_path = None
        self.cpu_num = None

        self.l_sample = []
        self.l_bam_path = []
        self.l_gtf_path = []
        self.l_tsv_path = []

        self.columns = []
        self.data = []
        self.all_refseq_id = None
        self.mg = mygene.MyGeneInfo()
        self.d_search_mygene = defaultdict(lambda: defaultdict(lambda: None))
        self.d_search_entrez = defaultdict(lambda: defaultdict(lambda: None))

    def start(self):
        self.get_args()
        try:
            self.check_file_path()
            self.echo_input_params()
            self._log("=== Analysis start. ===")
            self.run_stringtie()
            self.format_gtf_to_tsv()
            self.merge_tsv()
            self.get_entrez_info()
            self.get_go_info()
            self.dump_tsv()
            self._log("=== Analysis finish. ===")
        except:
            traceback.print_exc()
            sys.exit(1)

        return True

    def get_args(self):
        argparse_description = \
            "Expression analysis pipeline created by suecharo."
        help_bam = "Input all sorted bam file."
        help_gtf = "Input gtf file of hg19 or mm9."
        help_output = "Enter output dir.(default=./output)"
        help_cpu = "Input cpu num.(default=4)"

        parser = argparse.ArgumentParser(description=argparse_description)
        parser.add_argument("bam", nargs="+", help=help_bam)
        parser.add_argument("-g", nargs=1, required=True, metavar="GTF",
                            help=help_gtf)
        parser.add_argument("-o", nargs=1, metavar="OUTPUT", help=help_output)
        parser.add_argument("-p", nargs=1, type=int, metavar="CPU",
                            help=help_cpu)
        args = parser.parse_args()

        self.l_bam_path = args.bam
        self.gtf_path = args.g[0]
        self.output_dir_path = args.o
        if self.output_dir_path is None:
            self.output_dir_path = "./output"
        else:
            self.output_dir_path = self.output_dir_path[0]
        self.cpu_num = args.p
        if self.cpu_num is None:
            self.cpu_num = 4
        else:
            self.cpu_num = self.cpu_num[0]

        return True

    def check_file_path(self):
        l_tmp = []
        for bam_path in self.l_bam_path:
            if os.path.exists(bam_path):
                bam_path = os.path.abspath(bam_path)
                filename = os.path.basename(bam_path)
                l_filename = filename.split(".")
                sample_name = l_filename[0]
                ext = l_filename[-1]
                if ext == "bam":
                    pass
                else:
                    msg = "{} is wrong format.".format(bam_path)
                    raise ValueError(msg)
                l_tmp.append(bam_path)
                self.l_sample.append(sample_name)
            else:
                msg = "{} is not found.".format(bam_path)
                raise FileNotFoundError(msg)
        self.l_bam_path = l_tmp

        if os.path.exists(self.gtf_path):
            self.gtf_path = os.path.abspath(self.gtf_path)
            filename = os.path.basename(self.gtf_path)
            l_filename = filename.split(".")
            ext = l_filename[-1]
            if ext == "gtf":
                pass
            else:
                msg = "{} is wrong format.".format(self.gtf_path)
                raise ValueError(msg)
        else:
            raise FileNotFoundError("{} is not found.".format(self.gtf_path))

        if os.path.exists(self.output_dir_path):
            self.output_dir_path = os.path.abspath(self.output_dir_path)
        else:
            os.makedirs(os.path.abspath(self.output_dir_path))
            self.output_dir_path = os.path.abspath(self.output_dir_path)

        return True

    def echo_input_params(self):
        self._log("=== Input Param ===")
        msg = "Bam file : {}".format(" ".join(self.l_bam_path))
        self._log(msg)
        self._log("GTF file : {}".format(self.gtf_path))
        self._log("Output Dir : {}".format(self.output_dir_path))
        self._log("CPU Num : {}".format(self.cpu_num))

        return True


if __name__ == "__main__":
    my_expression_analysis = ExpressionAnalysisBam()
    my_expression_analysis.start()

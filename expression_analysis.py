#!/bin/env python
# coding: utf-8
import argparse
import os
import shlex
import subprocess
import sys
import traceback
from collections import defaultdict
from datetime import datetime

import yaml

import mygene


class ExpressionAnalysis(object):
    CONF_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__),
                             "conf.yml"))

    def __init__(self):
        self.trimmomatic_path = None
        self.adp_pe_path = None
        self.adp_se_path = None

        self.l_fastq = []
        self.fasta_path = None
        self.hisat2_index_path = None
        self.gtf_path = None
        self.output_dir_path = None
        self.cpu_num = None

        self.l_fastq_single = []
        self.l_fastq_paired = []
        self.l_sample_single = []
        self.l_sample_paired = []

        self.l_sample = []
        self.l_sam_path = []
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
            self.get_entrez_info()
            self.get_go_info()
            self.dump_tsv()
            self._log("=== Analysis finish. ===")
        except:
            traceback.print_exc()
            sys.exit(1)

        return True

    def _log(self, msg):
        now = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        print("[{}] {}".format(now, msg), flush=True)

        return True

    def get_args(self):
        argparse_description = \
            "Expression analysis pipeline created by suecharo."
        help_fastq = "Input all fastq file.(.fastq or .fa or .fa.gz) " + \
                     "(single-end=SRR000001.fq, " + \
                     "paired-end=SRR000002_1.fq,SRR000002_2.fq)"
        help_fasta = "Input fasta file of .fa or .fasta."
        help_index = "Input HISAT2 index path.(e.g. ./hisat2-index/hg19)"
        help_gtf = "Input gtf file of hg19 or mm9."
        help_output = "Enter output dir.(default=./output)"
        help_cpu = "Input cpu num.(default=4)"

        parser = argparse.ArgumentParser(description=argparse_description)
        parser.add_argument("fastq", nargs="+", help=help_fastq)
        reference_group = parser.add_mutually_exclusive_group()
        reference_group.add_argument("-a", nargs=1, metavar="FASTA",
                                     help=help_fasta)
        reference_group.add_argument("-i", nargs=1, metavar="INDEX",
                                     help=help_index)
        parser.add_argument("-g", nargs=1, required=True, metavar="GTF",
                            help=help_gtf)
        parser.add_argument("-o", nargs=1, metavar="OUTPUT", help=help_output)
        parser.add_argument("-p", nargs=1, type=int, metavar="CPU",
                            help=help_cpu)
        args = parser.parse_args()

        self.l_fastq = args.fastq
        self.fasta_path = args.a
        if self.fasta_path is not None:
            self.fasta_path = self.fasta_path[0]
        self.hisat2_index_path = args.i
        if self.hisat2_index_path is not None:
            self.hisat2_index_path = self.hisat2_index_path[0]
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

    def check_fastq_path(self):
        for fastq in self.l_fastq:
            if fastq.count(",") == 0:   # single end
                if os.path.exists(fastq) is False:
                    raise FileNotFoundError("{} is not found.".format(fastq))
                fastq = os.path.abspath(fastq)
                filename = os.path.basename(fastq)
                l_filename = filename.split(".")
                if len(l_filename) == 2:
                    base, ext = l_filename
                    if ext == "fq" or ext == "fastq":
                        pass
                    else:
                        raise ValueError("{} is wrong format.".format(fastq))
                elif len(l_filename) >= 3:
                    base = l_filename[0]
                    ext_1, ext_2 = l_filename[-2], l_filename[-1]
                    if (ext_1 == "fq" and ext_2 == "gz") or \
                            (ext_1 == "fastq" and ext_2 == "gz"):
                        pass
                    else:
                        raise ValueError("{} is wrong format.".format(fastq))
                else:
                    raise ValueError("{} is wrong format.".format(fastq))
                self.l_fastq_single.append(fastq)
                self.l_sample_single.append(base)

            elif fastq.count(",") == 1:   # paired end
                l_fastq = []
                l_base = []
                for fastq in fastq.split(","):
                    if os.path.exists(fastq) is False:
                        msg = "{} is not found.".format(fastq)
                        raise FileNotFoundError(msg)
                    fastq = os.path.abspath(fastq)
                    l_fastq.append(fastq)
                    filename = os.path.basename(fastq)
                    l_filename = filename.split(".")
                    if len(l_filename) == 2:
                        base, ext = l_filename
                        if ext == "fq" or ext == "fastq":
                            pass
                        else:
                            msg = "{} is not found.".format(fastq)
                            raise ValueError(msg)
                    elif len(l_filename) >= 3:
                        base = l_filename[0]
                        ext_1, ext_2 = l_filename[-2], l_filename[-1]
                        if (ext_1 == "fq" and ext_2 == "gz") or \
                                (ext_1 == "fastq" and ext_2 == "gz"):
                            pass
                        else:
                            msg = "{} is not found.".format(fastq)
                            raise ValueError(msg)
                    else:
                        raise ValueError("{} is wrong format.".format(fastq))
                    l_base.append(base)
                self.l_fastq_paired.append(l_fastq)
                self.l_sample_paired.append(l_base)

            else:
                raise FileNotFoundError("{} is not found.".format(fastq))

        return True

    def check_file_path(self):
        if self.fasta_path is not None:
            if os.path.exists(self.fasta_path):
                self.fasta_path = os.path.abspath(self.fasta_path)
                filename = os.path.basename(self.fasta_path)
                l_filename = filename.split(".")
                ext = l_filename[-1]
                if ext == "fa" or ext == "fasta":
                    pass
                else:
                    msg = "{} is wrong format.".format(self.fasta_path)
                    raise ValueError(msg)
            else:
                msg = "{} is not found.".format(self.fasta_path)
                raise FileNotFoundError(msg)

        if self.hisat2_index_path is not None:
            self.hisat2_index_path = os.path.abspath(self.hisat2_index_path)
            dirname = os.path.dirname(self.hisat2_index_path)
            if os.path.exists(dirname):
                pass
            else:
                msg = "{} is not found.".format(self.hisat2_index_path)
                raise FileNotFoundError(msg)

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
        if len(self.l_fastq_single) != 0:
            msg = "Single-end Fastq : {}".format(" ".join(self.l_fastq_single))
            self._log(msg)
        if len(self.l_fastq_paired) != 0:
            l_tmp = []
            for paired in self.l_fastq_paired:
                l_tmp.append(" ".join(paired))
            msg = "Paired-end Fastq : {}".format(" ".join(l_tmp))
            self._log(msg)
        if self.fasta_path is not None:
            self._log("Fasta : {}".format(self.fasta_path))
        if self.hisat2_index_path is not None:
            self._log("HISAT2 Index : {}".format(self.hisat2_index_path))
        self._log("GTF file : {}".format(self.gtf_path))
        self._log("Output Dir : {}".format(self.output_dir_path))
        self._log("CPU Num : {}".format(self.cpu_num))

        return True

    def read_trimmomatic_path(self):
        with open(ExpressionAnalysis.CONF_PATH, "r") as f:
            data = yaml.load(f)
            self.trimmomatic_path = data["TRIMMOMATIC_PATH"]
            self.adp_pe_path = data["TRIM_ADAPTER_PE"]
            self.adp_se_path = data["TRIM_ADAPTER_SE"]

        return True

    def _cmd_wrapper(self, cmd):
        self._log(cmd)
        l_cmd = shlex.split(cmd)
        proc = subprocess.Popen(l_cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        status = proc.returncode
        stdout = stdout.decode("utf-8").splitlines()
        stderr = stderr.decode("utf-8").splitlines()
        if status == 0:
            for message in stdout:
                self._log(message)
        else:
            for message in stderr:
                self._log(message)
                raise RuntimeError("{} is error.".format(cmd))

        return True

    def run_hisat2_build(self):
        self._log("Hisat2 build start.")
        output_dir = os.path.join(self.output_dir_path, "HISAT2_index")
        os.mkdir(output_dir)
        fastq_base = os.path.basename(self.fasta_path)
        fastq_base = ".".join(fastq_base.split(".")[:-1])
        output_path = os.path.join(output_dir, fastq_base)
        cmd = "hisat2-build {} {}".format(self.fasta_path, output_path)
        self._cmd_wrapper(cmd)
        self.hisat2_index_path = output_path

        return True

    def run_trimmomatic(self):
        self._log("Trimmomatic start.")
        output_dir = os.path.join(self.output_dir_path, "Trimmomatic")
        os.mkdir(output_dir)
        if len(self.l_sample_single) != 0:
            for i in range(len(self.l_sample_single)):
                sample = self.l_sample_single[i]
                fastq_path = self.l_fastq_single[i]
                length = self._check_fastq_length(fastq_path)
                output_path = \
                    os.path.join(output_dir, "{}_trimed.fastq".format(sample))
                l_cmd = ["java",
                         "-jar",
                         self.trimmomatic_path,
                         "SE",
                         "-threads",
                         self.cpu_num,
                         fastq_path,
                         output_path,
                         "ILLUMINACLIP:{}:2:40:15".format(self.adp_se_path),
                         "LEADING:20",
                         "TRAILING:20",
                         "SLIDINGWINDOW:4:15",
                         "MINLEN:{}".format(int(round(length * 0.8)))
                         ]
                cmd = " ".join(map(str, l_cmd))
                self._cmd_wrapper(cmd)
                self.l_fastq_single[i] = output_path

        if len(self.l_sample_paired) != 0:
            for i in range(len(self.l_sample_paired)):
                sample_1, sample_2 = self.l_sample_paired[i]
                fastq_path_1, fastq_path_2 = self.l_fastq_paired[i]
                length = self._check_fastq_length(fastq_path_1)
                output_1 = os.path.join(output_dir,
                                        "{}_trimed_1.fastq".format(sample_1))
                output_2 = os.path.join(output_dir,
                                        "{}_trimed_2.fastq".format(sample_2))
                tmp_1 = os.path.join(output_dir,
                                     "{}_tmp_1.fastq".format(sample_1))
                tmp_2 = os.path.join(output_dir,
                                     "{}_tmp_2.fastq".format(sample_2))
                l_cmd = ["java",
                         "-jar",
                         self.trimmomatic_path,
                         "PE",
                         "-threads",
                         self.cpu_num,
                         fastq_path_1,
                         fastq_path_2,
                         output_1,
                         tmp_1,
                         output_2,
                         tmp_2,
                         "ILLUMINACLIP:{}:2:40:15".format(self.adp_pe_path),
                         "LEADING:20",
                         "TRAILING:20",
                         "SLIDINGWINDOW:4:15",
                         "MINLEN:{}".format(length)
                         ]
                cmd = " ".join(map(str, l_cmd))
                self._cmd_wrapper(cmd)
                self.l_fastq_paired[i][0] = output_1
                self.l_fastq_paired[i][1] = output_2
                os.remove(tmp_1)
                os.remove(tmp_2)

        return True

    def _check_fastq_length(self, fastq_path):
        filename = os.path.basename(fastq_path)
        ext = filename.split(".")[-1]
        if ext == "gz":
            cmd = "gunzip -c {}".format(fastq_path)
            l_cmd = shlex.split(cmd)
            tmp_path = os.path.join(self.output_dir_path, "tmp.fq")
            with open(tmp_path, "w") as f:
                proc = subprocess.Popen(l_cmd, stdout=f)
            proc.wait()
            with open(tmp_path, "r") as f:
                for i, line in enumerate(iter(f.readline, "")):
                    if i == 1:
                        length = len(line) - 1
                        break
            os.remove(tmp_path)
        else:
            with open(fastq_path, "r") as f:
                for i, line in enumerate(iter(f.readline, "")):
                    if i == 1:
                        length = len(line) - 1
                        break

        return length

    def run_hisat2(self):
        self._log("HISAT2 start.")
        output_dir = os.path.join(self.output_dir_path, "HISAT2_sam")
        os.mkdir(output_dir)
        if len(self.l_sample_single) != 0:
            for i in range(len(self.l_sample_single)):
                fastq_path = self.l_fastq_single[i]
                sample = self.l_sample_single[i]
                output_path = os.path.join(output_dir, "{}.sam".format(sample))
                l_cmd = ["hisat2",
                         "-p",
                         self.cpu_num,
                         "--dta",
                         "-x",
                         self.hisat2_index_path,
                         "-U",
                         fastq_path,
                         "-S",
                         output_path
                         ]
                cmd = " ".join(map(str, l_cmd))
                self._cmd_wrapper(cmd)
                self.l_sample.append(sample)
                self.l_sam_path.append(output_path)

        if len(self.l_sample_paired) != 0:
            for i in range(len(self.l_sample_paired)):
                sample_1, sample_2 = self.l_sample_paired[i]
                sample = ""
                for j in range(min(len(sample_1), len(sample_2))):
                    if sample_1[j] == sample_2[j]:
                        sample += sample_1[j]
                    else:
                        break
                fastq_path_1, fastq_path_2 = self.l_fastq_paired[i]
                output_path = os.path.join(output_dir, "{}.sam".format(sample))
                l_cmd = ["hisat2",
                         "-p",
                         self.cpu_num,
                         "--dta",
                         "-x",
                         self.hisat2_index_path,
                         "-1",
                         fastq_path_1,
                         "-2",
                         fastq_path_2,
                         "-S",
                         output_path
                         ]
                cmd = " ".join(map(str, l_cmd))
                self._cmd_wrapper(cmd)
                self.l_sample.append(sample)
                self.l_sam_path.append(output_path)

        l_sample_sorted = sorted(self.l_sample)
        l_sam_path_sorted = []
        for sample in l_sample_sorted:
            ind = self.l_sample.index(sample)
            l_sam_path_sorted.append(self.l_sam_path[ind])
        self.l_sample = l_sample_sorted
        self.l_sam_path = l_sam_path_sorted

        return True

    def run_samtools_sort(self):
        self._log("Samtools sort start.")
        output_dir = os.path.join(self.output_dir_path, "HISAT2_bam")
        os.mkdir(output_dir)
        for i in range(len(self.l_sample)):
            sample = self.l_sample[i]
            sam_path = self.l_sam_path[i]
            output_path = os.path.join(output_dir, "{}.bam".format(sample))
            l_cmd = ["samtools",
                     "sort",
                     "-@",
                     self.cpu_num,
                     sam_path,
                     "-o",
                     output_path
                     ]
            cmd = " ".join(map(str, l_cmd))
            self._cmd_wrapper(cmd)
            self.l_bam_path.append(output_path)

        return True

    def run_stringtie(self):
        self._log("Stringtie start.")
        output_dir = os.path.join(self.output_dir_path, "Stringtie")
        os.mkdir(output_dir)
        for i in range(len(self.l_sample)):
            sample = self.l_sample[i]
            bam_path = self.l_bam_path[i]
            output_path = os.path.join(output_dir, "{}.gtf".format(sample))
            l_cmd = ["stringtie",
                     "-e",
                     bam_path,
                     "-G",
                     self.gtf_path,
                     "-o",
                     output_path,
                     "-p",
                     self.cpu_num
                     ]
            cmd = " ".join(map(str, l_cmd))
            self._cmd_wrapper(cmd)
            self.l_gtf_path.append(output_path)

        return True

    def format_gtf_to_tsv(self):
        self._log("Format GTF to tsv start.")
        output_dir = os.path.join(self.output_dir_path, "formated_tsv")
        os.mkdir(output_dir)
        l_columns = ["Chromosome", "Start", "End", "Width", "Strand",
                     "Gene_ID", "RefSeq_ID", "Coverage", "FPKM", "TPM"]
        for i in range(len(self.l_sample)):
            sample = self.l_sample[i]
            gtf_path = self.l_gtf_path[i]
            output_path = os.path.join(output_dir, "{}.tsv".format(sample))

            with open(output_path, "w") as f_w:
                all_write = []
                all_write.append("\t".join(l_columns))
                with open(gtf_path, "r") as f_r:
                    for line in iter(f_r.readline, ""):
                        if line[0] == "#" or line == "":
                            continue
                        l_line = line.split("\t")
                        if l_line[2] != "transcript":
                            continue
                        l_write = []
                        for i in [0, 3, 4]:
                            l_write.append(l_line[i])
                        l_write.append(str(int(l_line[4]) - int(l_line[3])))
                        l_write.append(l_line[6])
                        l_misc = l_line[8].split(";")
                        for i in [0, 1, 3, 4, 5]:
                            l_write.append(l_misc[i].split('"')[1])
                        all_write.append("\t".join(l_write))

                f_w.write("\n".join(all_write))
            self.l_tsv_path.append(output_path)

        return True

    def merge_tsv(self):
        self._log("Merge tsv start.")
        l_dict = []
        tmp_columns = ["Chromosome", "Start", "End", "Width", "Strand",
                       "Gene_ID", "Coverage", "FPKM", "TPM"]
        tmp_index = [0, 1, 2, 3, 4, 5, 7, 8, 9]
        d_all_tmp = defaultdict(lambda: defaultdict(lambda: None))

        for tsv_path in self.l_tsv_path:
            tmp_dict = defaultdict(lambda: None)
            with open(tsv_path, "r") as f:
                data = f.read()
                data = data.split("\n")
            for row in data[1:]:
                if row == "":
                    continue
                l_ele = row.split("\t")
                refseq_id = l_ele[6]
                tmp_dict[refseq_id] = defaultdict(lambda: None)
                for column, ind in zip(tmp_columns, tmp_index):
                    d_all_tmp[refseq_id][column] = l_ele[ind]
                    tmp_dict[refseq_id][column] = l_ele[ind]
            l_dict.append(tmp_dict)

        l_sort = []
        for key, value in d_all_tmp.items():
            l_sort.append([key, value["Chromosome"], value["Start"]])
        l_sort.sort(key=lambda x: (x[1], x[2]))
        self.all_refseq_id = [row[0] for row in l_sort]
        for refseq_id in self.all_refseq_id:
            self.data.append([refseq_id])

        self.columns.append("RefSeq_ID")
        for i in range(len(self.l_sample)):
            sample = self.l_sample[i]
            tmp_dict = l_dict[i]
            self.columns.append("{}_FPKM".format(sample))
            self.columns.append("{}_TPM".format(sample))
            self.columns.append("{}_Coverage".format(sample))
            for i, refseq_id in enumerate(self.all_refseq_id):
                if tmp_dict[refseq_id] is None:
                    self.data[i].append(0)
                    self.data[i].append(0)
                    self.data[i].append(0)
                else:
                    self.data[i].append(tmp_dict[refseq_id]["FPKM"])
                    self.data[i].append(tmp_dict[refseq_id]["TPM"])
                    self.data[i].append(tmp_dict[refseq_id]["Coverage"])

        l_columns = ["Gene_ID", "Chromosome", "Start", "End", "Width",
                     "Strand"]
        for column in l_columns:
            self.columns.append(column)
            if column in ["Start", "End", "Width"]:
                for i, refseq_id in enumerate(self.all_refseq_id):
                    self.data[i].append(int(d_all_tmp[refseq_id][column]))
            else:
                for i, refseq_id in enumerate(self.all_refseq_id):
                    self.data[i].append(d_all_tmp[refseq_id][column])

        return True

    def get_entrez_info(self):
        self._log("Get entrez info start.")
        len_l_all_refseq = len(self.all_refseq_id)
        self._log("Length list of RefSeq_ID is {}".format(len_l_all_refseq))
        self.columns.append("Entrez_Gene_ID")
        self.columns.append("Description")
        ind_gene_id = self.columns.index("Gene_ID")
        for i in range(len_l_all_refseq):
            row = self.data[i]
            refseq_id = self.all_refseq_id[i]
            gene_id = row[ind_gene_id]
            self.d_search_mygene[refseq_id]["Gene_ID"] = gene_id
        self._querymany_refseq_id()

        for i in range(len_l_all_refseq):
            row = self.data[i]
            refseq_id = self.all_refseq_id[i]
            entrez_gene_id = self.d_search_mygene[refseq_id]["Entrez_Gene_ID"]
            description = self.d_search_mygene[refseq_id]["Description"]
            self.data[i].append(entrez_gene_id)
            self.data[i].append(description)

        return True

    def _querymany_refseq_id(self):
        d_ret = self.mg.querymany(self.all_refseq_id, scopes="refseq",
                                  returnall=True)
        l_out = d_ret["out"]
        l_missing = d_ret["missing"]
        for out in l_out:
            refseq_id = out["query"]
            if "entrezgene" in out and "name" in out and "symbol" in out:
                entrez_gene_id = out["entrezgene"]
                description = out["name"]
                gene_id = out["symbol"]
                val_gene_id = self.d_search_mygene[refseq_id]["Gene_ID"]
                if gene_id == val_gene_id:
                    self.d_search_mygene[refseq_id]["Entrez_Gene_ID"] = \
                        entrez_gene_id
                    self.d_search_mygene[refseq_id]["Description"] = \
                        description
                else:
                    l_missing.append(refseq_id)
            else:
                l_missing.append(refseq_id)

        d_refind = dict()
        l_nothing = []
        for refseq_id in l_missing:
            l_refseq_id = refseq_id.split("_")
            if len(l_refseq_id) <= 2:
                l_nothing.append(refseq_id)
            else:
                short_refseq_id = "_".join(l_refseq_id[:2])
                d_refind[short_refseq_id] = refseq_id

        d_ret = self.mg.querymany(list(d_refind.keys()), scopes="refseq",
                                  returnall=True)
        l_out = d_ret["out"]
        l_missing = d_ret["missing"]
        for short_refseq_id in l_missing:
            l_nothing.append(d_refind[short_refseq_id])
        for out in l_out:
            short_refseq_id = out["query"]
            refseq_id = d_refind[short_refseq_id]
            if "entrezgene" in out and "name" in out and "symbol" in out:
                entrez_gene_id = out["entrezgene"]
                description = out["name"]
                gene_id = out["symbol"]
                val_gene_id = self.d_search_mygene[refseq_id]["Gene_ID"]
                if gene_id == val_gene_id:
                    self.d_search_mygene[refseq_id]["Entrez_Gene_ID"] = \
                        entrez_gene_id
                    self.d_search_mygene[refseq_id]["Description"] = \
                        description
                else:
                    l_nothing.append(refseq_id)
            else:
                l_nothing.append(refseq_id)

        d_nothing = dict()
        for refseq_id in l_nothing:
            gene_id = self.d_search_mygene[refseq_id]["Gene_ID"]
            d_nothing[gene_id] = refseq_id

        d_ret = self.mg.querymany(list(d_nothing.keys()), scopes="symbol",
                                  size=1, returnall=True)
        l_out = d_ret["out"]
        s_missing = set(d_ret["missing"])
        for out in l_out:
            gene_id = out["query"]
            refseq_id = d_nothing[gene_id]
            if "entrezgene" in out and "name" in out:
                entrez_gene_id = out["entrezgene"]
                description = out["name"]
                self.d_search_mygene[refseq_id]["Entrez_Gene_ID"] = \
                    entrez_gene_id
                self.d_search_mygene[refseq_id]["Description"] = \
                    description
            else:
                if "notfound" not in out:
                    s_missing.add(gene_id)

        for i, gene_id in enumerate(list(s_missing)):
            self._log("Query each refseq id {} / {}".format(i + 1,
                                                            len(s_missing)))
            ret = self.mg.query(gene_id, size=1)
            refseq_id = d_nothing[gene_id]
            hits = ret["hits"]
            if len(hits) == 1:
                ele = hits[0]
                if "name" in ele:
                    description = ele["name"]
                    self.d_search_mygene[refseq_id]["Description"] = \
                        description
                if "entrezgene" in ele:
                    entrez_gene_id = ele["entrezgene"]
                    self.d_search_mygene[refseq_id]["Entrez_Gene_ID"] = \
                        entrez_gene_id

        return True

    def get_go_info(self):
        self._log("Get GO info start.")
        self.columns.append("Cellular_Component_Accession")
        self.columns.append("Cellular_Component_Name")
        self.columns.append("Molecular_Function_Accession")
        self.columns.append("Molecular_Function_Name")
        self.columns.append("Biological_Process_Accession")
        self.columns.append("Biological_Process_Name")

        self._querymany_entrez_gene_id()
        ind_entrez_gene_id = self.columns.index("Entrez_Gene_ID")
        for i in range(len(self.all_refseq_id)):
            row = self.data[i]
            entrez_gene_id = str(row[ind_entrez_gene_id])
            for category in ["CC", "MF", "BP"]:
                for key in ["_accession", "_term"]:
                    d_key = category + key
                    ele = self.d_search_entrez[entrez_gene_id][d_key]
                    self.data[i].append(ele)

        return True

    def _querymany_entrez_gene_id(self):
        s_entrez_gene_id = set()
        ind_entrez_gene_id = self.columns.index("Entrez_Gene_ID")
        for i in range(len(self.all_refseq_id)):
            row = self.data[i]
            entrez_gene_id = row[ind_entrez_gene_id]
            if entrez_gene_id is not None:
                s_entrez_gene_id.add(entrez_gene_id)

        l_ret = self.mg.getgenes(list(s_entrez_gene_id), fields="go")
        l_go_categories = ["CC", "MF", "BP"]
        for ret in l_ret:
            entrez_gene_id = str(ret["query"])
            d_tmp = defaultdict(list)
            if "go" in ret:
                d_go = ret["go"]
                for category in l_go_categories:
                    if category in d_go:
                        l_go = d_go[category]
                        for ele_go in l_go:
                            if isinstance(ele_go, dict):
                                if "id" in ele_go:
                                    accession = ele_go["id"]
                                    key = category + "_accession"
                                    d_tmp[key].append(accession)
                                if "term" in ele_go:
                                    name = ele_go["term"]
                                    key = category + "_term"
                                    d_tmp[key].append(name)
            for category in l_go_categories:
                for key in ["_accession", "_term"]:
                    d_key = category + key
                    l_ele = d_tmp[d_key]
                    if len(l_ele) == 0:
                        ele = None
                    else:
                        ele = "//".join(l_ele)
                    self.d_search_entrez[entrez_gene_id][d_key] = ele

        return True

    def dump_tsv(self):
        self._log("Start dump tsv.")
        output_path = os.path.join(self.output_dir_path, "analysed.tsv")
        l_write = []
        l_write.append("\t".join(map(str, self.columns)))
        for row in self.data:
            l_write.append("\t".join(map(str, row)))
        s_write = "\n".join(l_write)
        with open(output_path, "w") as f:
            f.write(s_write)

        return True


if __name__ == "__main__":
    my_expression_analysis = ExpressionAnalysis()
    my_expression_analysis.start()

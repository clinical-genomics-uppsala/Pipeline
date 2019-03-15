# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from enum import Enum, auto, unique
import enum
import re


_cds_pattern =  re.compile(r'^c\..+|^-$')
_aa_pattern =  re.compile(r'^p\..+|^-$')
_exon_intron_pattern =  re.compile(r'^exon\d+$|^intronic$')
_chr_pattern =  re.compile(r'^chr\d+$|^\d+$')
_nc_pattern =  re.compile(r'^NC_0+\d+\.\d+$')


@unique
class VariantClass(Enum):
    hotspot = 1
    indel  =  2
    check =  3
    other = 4

    def __str__(self):
        return '%s-%s' % (self.value, self.name)

@unique
class ReportClass(Enum):
    hotspot = auto()
    indel = auto()
    region = auto()
    region_all = auto()

    def __str__(self):
        return '%s' % self.name


class PileupPosition(object):
    def __init__(self, CHROMOSOME, POSITION, REFERENCE_BASE, DEPTH):
        self.CHROMOSOME = CHROMOSOME
        self.POSITION = int(POSITION)
        self.REFERENCE_BASE = REFERENCE_BASE
        self.DEPTH = int(DEPTH)


class Hotspot(object):
    def __init__(self, CHROMOSOME, START, END, GENE, CDS_MUTATION_SYNTAX, AA_MUTATION_SYNTAX, REPORT, COMMENT, EXON, ACCESSION_NUMBER):
        self.CHROMOSOME = CHROMOSOME
        self.START = START
        self.END = END
        self.GENE = GENE
        self.CDS_MUTATION_SYNTAX = CDS_MUTATION_SYNTAX
        self.AA_MUTATION_SYNTAX = AA_MUTATION_SYNTAX
        self.REPORT = REPORT
        self.COMMENT = COMMENT
        self.EXON = EXON
        self.ACCESSION_NUMBER = ACCESSION_NUMBER

        if not isinstance(self.START, int):
            raise ValueError("Start position should be an integer: %s!" % self.START)

        if not isinstance(self.END, int):
            raise ValueError("End position should be an integer: %s!" % self.END)

        if self.START > self.END:
            raise ValueError("Start cordinte cannot be larget then stop coordinate! % > %" % (self.START, self.END))

        if not _cds_pattern.match(self.CDS_MUTATION_SYNTAX):
            raise ValueError("Incorrect cds syntax %s! Should start with \"c.\" or set to \"-\" if empty." % self.CDS_MUTATION_SYNTAX)

        if not _aa_pattern.match(self.AA_MUTATION_SYNTAX):
            raise ValueError("Incorrect aa syntax: %s! Should start with \"p.\" or set to \"-\" if empty." % self.AA_MUTATION_SYNTAX)

        if self.REPORT in ReportClass.__members__:
            raise ValueError("report value (%s) not found in  Enum class %s!" % (self.REPORT,list(ReportClass)))

        if not _exon_intron_pattern.match(self.EXON):
            raise ValueError("Exon value should have the following format: exon or intronic. not %" % self.EXON)

        self.TOTAL_DEPTH = ["-"] * (self.END - self.START + 1)
        self.VARIANTS = []


    def check_overlapp(self, start,stop=None):
        return  (stop and self.START <= stop <= self.END) or (self.START <= start <= self.END)


    def add_variant(self, variant):
        self.variants.append(variant)

class ChrTranslater(object):
    def __init__(self, mapper_file):
        self.chr_to_nc = dict()
        self.nc_to_chr = dict()

        with open(mapper_file) as mapping:
            for line in mapping:
                if not line.startswith("#"):
                    columns = line.split("\t")
                    if _chr_pattern.match(columns[0]) and _nc_pattern.match(columns[1]):
                        self.chr_to_nc[columns[0]] = columns[1]
                        self.nc_to_chr[columns[1]] = columns[0]
                    elif _chr_pattern.match(columns[1]) and _nc_pattern.match(columns[0]):
                        self.chr_to_nc[columns[1]] = columns[0]
                        self.nc_to_chr[columns[0]] = columns[1]
                    else:
                        raise Exception("Unexpected column values for column 1/2")

    def get_chr_value(self, nc_id):
        return self.nc_to_chr[nc_id]

    def get_nc_value(self, chr_id):
        return self.chr_to_nc[chr_id]

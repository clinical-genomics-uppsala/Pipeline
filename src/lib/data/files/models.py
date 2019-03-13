# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from enum import Enum, auto, unique
import enum
import re

import pysam

import src.lib.data.files.vcf as vcf

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

    def __str__(self):
        return "{}\t{}\t{}\t{}".format(self.CHROMOSOME, self.POSITION, self.REFERENCE_BASE, self.DEPTH)


class MultiBpVariant(object):
    def __init__(self,CHROMOSOME, START, STOP, REFERENCE, VARIANT, GENE, CDS_CHANGE, AA_CHANGE, TRANSCRIPT):
        self.CHROMOSOME = CHROMOSOME
        self.START = START
        self.STOP = STOP
        self.REFERENCE = REFERENCE
        self.VARIANT = VARIANT
        self.GENE = GENE
        self.CDS_CHANGE = CDS_CHANGE
        self.AA_CHANGE = AA_CHANGE
        self.TRANSCRIPT = TRANSCRIPT.split(".")[0]

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}".format(self.CHROMOSOME, self.START, self.STOP, self.REFERENCE, self.VARIANT)

class MultiBpVariantData(object):
    def __init__(self, data_file):
        self.data = {}
        from src.lib.data.files.multibp import Reader
        reader = Reader(data_file)
        for bp in reader:
            key = "{}:{}:{}:{}:{}".format(bp.CHROMOSOME, bp.START, bp.STOP, bp.REFERENCE, bp.VARIANT)
            if key in self.data:
                raise Exception("Trying to overwrite data för: " + key)
            self.data[key] = bp


    def get_data(self, chromosome, start, stop, reference, variant):
        return self.data.get("{}:{}:{}:{}:{}".format(chromosome, start, stop, reference, variant), None)

class Hotspot(object):
    def __init__(self, CHROMOSOME, START, END, GENE, CDS_MUTATION_SYNTAX, AA_MUTATION_SYNTAX, REPORT, COMMENT, EXON, ACCESSION_NUMBER, PRINT_ALL=False):
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

        self.DEPTH_VARIANTS = [{'depth': "-", 'extended': False, 'variants': []} for i in range((self.END - self.START + 1))]

        self.EXTENDED_START = self.START
        self.EXTENDED_END = self.END
        self.PRINT_ALL = PRINT_ALL
        self.VARIANT_ADDED = False

    def check_overlapp(self, chrom, region_start, region_stop, start,stop=None):
        #print(self.CHROMOSOME + " == " + chrom + " and (( " + str(stop) + " is not None and " + str(region_start) + " <= " + str(stop) + " and " + str(start) + " <= " + str(region_stop) + ") or (" + str(region_start) + " <= " + str(start) + " <= " + str(region_stop) + ")) " + str(self.CHROMOSOME == chrom) + " " + str((stop is not None and region_start <= stop and start <= region_stop)))
        return  self.CHROMOSOME == chrom and ((stop is not None and region_start <= stop and start <= region_stop) or (region_start <= start <= region_stop))


    def add_variant(self, variant, chr_translater):
        if isinstance(variant, pysam.VariantRecord):
            v_start = variant.start + 1
            v_stop = variant.stop + 1
            #print(str(variant))
            if self.check_overlapp(chr_translater.get_nc_value(variant.chrom), self.START, self.END, v_start , v_stop):
                if self.EXTENDED_END < v_stop or v_start < self.EXTENDED_START:
                    if variant.start < self.EXTENDED_START or self.EXTENDED_END < variant.stop:
                        new_start = min(v_start,self.EXTENDED_START)
                        new_end = max(v_stop, self.EXTENDED_END)
                        new_depth_var = []
                        for i in range(new_end - new_start + 1):
                            if self.START <= i + new_start <= self.END:
                                new_depth_var.append(self.DEPTH_VARIANTS[new_start - self.EXTENDED_START + i])
                            else:
                                new_depth_var.append({'depth': "-", 'extended': True, 'variants': []})
                        self.DEPTH_VARIANTS = new_depth_var
                        self.EXTENDED_START = new_start
                        self.EXTENDED_END = new_end
                position = self.START - self.EXTENDED_START
                try:
                    self.DEPTH_VARIANTS[position]['variants'].append(variant)
                except:
                    self.DEPTH_VARIANTS[position]['variants'] = [variant]
                self.VARIANT_ADDED = True
                return True
        return False

    def add_depth(self, depth, chr_translater):
        if isinstance(depth, PileupPosition):
            #print("Checking: " + str(depth)+ "\t" + str(depth.POSITION))
            if self.check_overlapp(chr_translater.get_nc_value(depth.CHROMOSOME), self.EXTENDED_START, self.EXTENDED_END, depth.POSITION, None):
                #print("Adding depth: " + str(depth)+ "\t" + str(depth.POSITION))
                #print(depth.POSITION-self.EXTENDED_START)
                self.DEPTH_VARIANTS[depth.POSITION-self.EXTENDED_START]['depth'] = int(depth.DEPTH)
                return True
        return False

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

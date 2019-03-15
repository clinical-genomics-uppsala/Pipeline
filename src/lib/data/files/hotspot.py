# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import codecs
from collections import OrderedDict
import gzip
import re

from .models import Hotspot, ReportClass

class Reader(object):

    def __init__(self, filename, compressed=None, strict_whitespace=False, encoding='ascii'):
        super().__init__()

        if not filename:
            raise Exception("A filepath needs to be specified: %" % filepath)

        if compressed:
            self._reader = codecs.getreader(encoding)(self._reader)
        else:
            if compressed is None:
                compressed = filename.endswith('.gz')
            self._reader = open(filename, 'rb' if compressed else 'rt')
        self.filename = filename

        self._separator = '\t'

        self._row_pattern =  re.compile(self._separator)

        self.reader = (line.strip() for line in self._reader if line.strip())

        self._parse_header()

        self.mapper = dict()


    def _parse_header(self):
        line = next(self.reader)
        if not line.startswith("#"):
            raise Exception("Missing header row!")
        line = line.lstrip("#")
        self.header = OrderedDict([(name.upper(),index) for index, name in enumerate(line.split("\t"))])


    def __iter__(self):
        return self


    def next(self):
        '''Return the next record in the file.'''
        line = next(self.reader)
        row = self._row_pattern.split(line.rstrip())
        chrom = row[self.header['CHR']]

        start = int(row[self.header['START']])
        stop = int(row[self.header['END']])
        gene = row[self.header['GENE']]

        cds = row[self.header['CDS_MUTATION_SYNTAX']]
        aa = row[self.header['AA_MUTATION_SYNTAX']]
        report = ReportClass[row[self.header['REPORT']]]
        comment = row[self.header['COMMENT']]
        exon = row[self.header['EXON']]
        accession_numbber = row[self.header['ACCESSION_NUMBER']]

        return Hotspot(chrom, start, stop, gene, cds, aa, report, comment, exon, accession_numbber)

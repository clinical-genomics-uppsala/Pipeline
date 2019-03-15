# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import re
import gzip

from .models import PileupPosition

class Reader(object):
    def __init__(self, filename, compressed=None, encoding='ascii'):
        super().__init__()

        if not filename:
            raise Exception("A filepath needs to be specified: %" % filepath)

        if compressed is None:
            compressed = filename.endswith('.gz')

        self._reader = gzip.open(filename, 'rt') if compressed else open(filename,'rt')

        self.filename = filename

        self.reader = (line.strip() for line in self._reader if line.strip())

        self._row_pattern =  re.compile("\t")

    def __del__(self):
        if self._reader and not self._reader.closed:
            self._reader.close()

    def __iter__(self):
        return self


    def next(self):
        return self.__next__()

    def __next__(self):
        '''Return the next record in the file.'''
        line = next(self.reader)
        row = self._row_pattern.split(line.rstrip())
        chrom = row[0]
        pos = row[1]
        ref_base = row[2]
        depth = row[3]

        return PileupPosition(chrom, pos, ref_base, depth)

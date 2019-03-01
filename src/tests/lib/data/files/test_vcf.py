# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import logging
import os
import shutil
import tempfile
import unittest

logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

class TestVcf(unittest.TestCase):
    def setUp(self):
        # create fixtures
        self.tempdir = tempfile.mkdtemp()
        with open(os.path.join(self.tempdir,"in.vcf"),'w') as hotspots:
            hotspots.write("##fileformat=VCFv4.0\n")
            hotspots.write("##contig=<ID=chr7,length=159138663,assembly=hg19>\n")
            hotspots.write("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n")
            hotspots.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            hotspots.write("chr7\t55249071\t.\tC\tT\t4399\tPASS\tDP4=17465,17450,117,118")


    def tearDown(self):
        # delete fixtures
        shutil.rmtree(self.tempdir)

    def test_add_AD_field_using(self):
        from src.lib.data.files.vcf import add_AD_field_using_DP4
        add_AD_field_using_DP4(os.path.join(self.tempdir,"in.vcf"), os.path.join(self.tempdir,"out.vcf"))

        with open(os.path.join(self.tempdir,"out.vcf"), 'r') as output_vcf:
            lines = output_vcf.readlines()
            self.assertEqual("fileformat" in lines[0], True)
            self.assertEqual("ID=PASS" in lines[1], True)
            self.assertEqual("ID=DP4" in lines[3], True)
            self.assertEqual("ID=AD" in lines[4], True)
            self.assertEqual("AD=34915,235" in lines[6], True)


if __name__ == '__main__':
    import logging
    import sys
    logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    unittest.main()

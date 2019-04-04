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
            hotspots.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
            hotspots.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            hotspots.write("chr7\t55249071\t.\tC\tT\t4399\tPASS\tDP=35150;DP4=17465,17450,117,118")
        with open(os.path.join(self.tempdir,"in.withoutcontigs.vcf"),'w') as hotspots:
            hotspots.write("##fileformat=VCFv4.0\n")
            hotspots.write("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n")
            hotspots.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
            hotspots.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allelic Frequency\">\n")
            hotspots.write("##INFO=<ID=AD,Number=.,Type=String,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n")
            hotspots.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            hotspots.write("chr7\t55249071\t.\tC\tT\t4399\tPASS\tDP=35150;DP4=17465,17450,117,118;AD=34915,235;AF=0.5")
        with open(os.path.join(self.tempdir,"reference_info"),'w') as hotspots:
            hotspots.write("#Chr name\tNC ID\tName\tLength\n")
            hotspots.write("chr7\tNC_000007.12\tChr7#NC_000007.12#1#158821424#-1\t158821424")


    def tearDown(self):
        # delete fixtures
        shutil.rmtree(self.tempdir)


    def test_add_AD_field_using(self):
        from src.lib.data.files.vcf import add_AD_field_using_DP4
        add_AD_field_using_DP4(os.path.join(self.tempdir,"in.vcf"), os.path.join(self.tempdir,"out.vcf"))
        with open(os.path.join(self.tempdir,"out.vcf"), 'r') as output_vcf:
            lines = output_vcf.readlines()
            self.assertTrue("fileformat" in lines[0])
            self.assertTrue("ID=PASS" in lines[1])
            self.assertTrue("contig" in lines[2])
            self.assertTrue("ID=DP4," in lines[3])
            self.assertTrue("ID=DP," in lines[4])
            self.assertTrue("ID=AD" in lines[5])
            self.assertTrue("AD=34915,235" in lines[7])

    def test_create_sample_format_from_info(self):
        from src.lib.data.files.vcf import create_sample_format_from_info_lofreq
        from src.lib.data.files.vcf import add_contigs_to_header
        add_contigs_to_header(
            os.path.join(self.tempdir,"in.withoutcontigs.vcf"),
            os.path.join(self.tempdir,"out.withcontigs.vcf"),
            os.path.join(self.tempdir,"reference_info"),
            "hg19")
        create_sample_format_from_info_lofreq(
            "sample1",
             os.path.join(self.tempdir,"out.withcontigs.vcf"),
             os.path.join(self.tempdir,"sample.vcf"))
        with open(os.path.join(self.tempdir,"sample.vcf"), 'r') as output_vcf:
             lines = output_vcf.readlines()
             self.assertTrue("fileformat" in lines[0])
             self.assertTrue("ID=PASS" in lines[1])
             self.assertTrue(lines[2].startswith("##INFO=<ID=DP4,"))
             self.assertTrue(lines[3].startswith("##INFO=<ID=DP,"))
             self.assertTrue(lines[4].startswith("##INFO=<ID=AF,"))
             self.assertTrue(lines[5].startswith("##INFO=<ID=AD,"))
             self.assertEqual(lines[6].rstrip(), "##contig=<ID=chr7,length=158821424,assembly=hg19>")
             self.assertTrue(lines[7].startswith("##FORMAT=<ID=AF,"))
             self.assertTrue(lines[8].startswith("##FORMAT=<ID=AD,"))
             self.assertTrue(lines[9].startswith("##FORMAT=<ID=DP,"))
             self.assertTrue(lines[10].startswith("##FORMAT=<ID=DP4,"))
             self.assertTrue(lines[11].startswith("##FORMAT=<ID=GT,"))
             self.assertEqual(lines[12], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
             self.assertEqual(lines[13].rstrip(), "chr7\t55249071\t.\tC\tT\t4399\tPASS\tDP=35150;DP4=17465,17450,117,118;AD=34915,235;AF=0.5\tGT:AF:DP4:DP:AD\t0/1:0.5:17465,17450,117,118:35150:34915,235")


    def test_add_contigs_to_header(self):
        from src.lib.data.files.vcf import add_contigs_to_header
        add_contigs_to_header(
            os.path.join(self.tempdir,"in.withoutcontigs.vcf"),
            os.path.join(self.tempdir,"out.withcontigs.vcf"),
            os.path.join(self.tempdir,"reference_info"),
            "hg19")

        with open(os.path.join(self.tempdir,"out.withcontigs.vcf"), 'r') as output_vcf:
            lines = output_vcf.readlines()
            self.assertTrue("fileformat" in lines[0])
            self.assertTrue("ID=PASS" in lines[1])
            self.assertTrue("INFO=<ID=DP4," in lines[2])
            self.assertTrue("INFO=<ID=DP," in lines[3])
            self.assertTrue("INFO=<ID=AF" in lines[4])
            self.assertTrue("INFO=<ID=AD" in lines[5])
            self.assertEqual(lines[6].rstrip(), "##contig=<ID=chr7,length=158821424,assembly=hg19>")



if __name__ == '__main__':
    import logging
    import sys
    logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    unittest.main()

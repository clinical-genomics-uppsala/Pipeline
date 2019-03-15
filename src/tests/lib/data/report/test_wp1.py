# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import logging
import os
import shutil
import tempfile
import unittest

logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

class TestWp1Reports(unittest.TestCase):
    def setUp(self):
        # create fixtures
        self.tempdir = tempfile.mkdtemp()
        with open(os.path.join(self.tempdir,"hotspot"),'w') as hotspots:
            hotspots.write("#Chr\tStart\tEnd\tGene\tCDS_mutation_syntax\tAA_mutation_syntax\tReport\tcomment\tExon\tAccession_number\n")
            hotspots.write("NC_000007.13\t55249003\t55249003\tEGFR\tc.C2301\tp.A767\thotspot\t-\texon20\tNM_005228\n")
            #hotspots.write("NC_000001.10\t115252202\t115252202\tNRAS\tc.C438\tp.A146\thotspot\t-\texon4\tNM_002524\n")
            #hotspots.write("NC_000017.10\t37880979\t37881164\tERBB2\t-\t-\tindel\t-\texon20\tNM_004448\n")
            #hotspots.write("NC_000012.11\t25398279\t25398279\tKRAS\tc.G40\tp.V14\tregion\t-\texon2\tNM_004985\n")
            #hotspots.write("NC_000017.10\t37880986\t37880987\tERBB2\t-\tp.Y772\tregion_all\t-\texon20\tNM_004448\n")

#M_005228	-	2-indel	yes	ok	4200	2996	1204	0.286666666667	-	-	-	No	ID=COSM12384;OCCURENCE=55(lung)	Tyrosine_kinase_inhibitor_response	drug-response	-	-	-	-	-	-	-	-	-	+1161|-43	-	-	NC_000007.13	55242467	55242485	AATTAAGAGAAGCAACATC	T	EGFR:NM_005228:exon19:c.2237_2255T
    def tearDown(self):
        # delete fixtures
        shutil.rmtree(self.tempdir)


    def test_header_parsing(self):
         from src.lib.data.files.hotspot import Reader
         reader = Reader(os.path.join(self.tempdir,"hotspot"))
         self.assertEqual(10, len(reader.header.keys()))
         self.assertEqual(0, reader.header["CHR"])
         self.assertEqual(1, reader.header["START"])
         self.assertEqual(2, reader.header["END"])
         self.assertEqual(3, reader.header["GENE"])
         self.assertEqual(4, reader.header["CDS_MUTATION_SYNTAX"])
         self.assertEqual(5, reader.header["AA_MUTATION_SYNTAX"])
         self.assertEqual(6, reader.header["REPORT"])
         self.assertEqual(7, reader.header["COMMENT"])
         self.assertEqual(8, reader.header["EXON"])
         self.assertEqual(9, reader.header["ACCESSION_NUMBER"])


    def test_record_parsing(self):
        from src.lib.data.files.hotspot import Reader
        from src.lib.data.files.models import ReportClass
        reader = Reader(os.path.join(self.tempdir,"hotspot"))

        record = reader.next()
        self.assertEqual("NC_000001.10", record.CHROMOSOME)
        self.assertEqual(115252190, record.START)
        self.assertEqual(115252349, record.END)
        self.assertEqual("NRAS", record.GENE)
        self.assertEqual("-", record.CDS_MUTATION_SYNTAX)
        self.assertEqual("-", record.AA_MUTATION_SYNTAX)
        self.assertEqual(ReportClass.indel, record.REPORT)
        self.assertEqual("-", record.COMMENT)
        self.assertEqual("exon4", record.EXON)
        self.assertEqual("NM_002524", record.ACCESSION_NUMBER)
        self.assertEqual(160, len(record.TOTAL_DEPTH))
        self.assertEqual(["-"]*160, record.TOTAL_DEPTH)

        record = reader.next()
        self.assertEqual("NC_000001.10", record.CHROMOSOME)
        self.assertEqual(115252202, record.START)
        self.assertEqual(115252202, record.END)
        self.assertEqual("NRAS", record.GENE)
        self.assertEqual("c.C438", record.CDS_MUTATION_SYNTAX)
        self.assertEqual("p.A146", record.AA_MUTATION_SYNTAX)
        self.assertEqual(ReportClass.hotspot, record.REPORT)
        self.assertEqual("-", record.COMMENT)
        self.assertEqual("exon4", record.EXON)
        self.assertEqual("NM_002524", record.ACCESSION_NUMBER)
        self.assertEqual(1, len(record.TOTAL_DEPTH))
        self.assertEqual(["-"], record.TOTAL_DEPTH)

        record = reader.next()
        self.assertEqual("NC_000017.10", record.CHROMOSOME)
        self.assertEqual(37880979, record.START)
        self.assertEqual(37881164, record.END)
        self.assertEqual("ERBB2", record.GENE)
        self.assertEqual("-", record.CDS_MUTATION_SYNTAX)
        self.assertEqual("-", record.AA_MUTATION_SYNTAX)
        self.assertEqual(ReportClass.indel, record.REPORT)
        self.assertEqual("-", record.COMMENT)
        self.assertEqual("exon20", record.EXON)
        self.assertEqual("NM_004448", record.ACCESSION_NUMBER)
        self.assertEqual(186, len(record.TOTAL_DEPTH))
        self.assertEqual(["-"]*186, record.TOTAL_DEPTH)

        record = reader.next()
        self.assertEqual("NC_000012.11", record.CHROMOSOME)
        self.assertEqual(25398279, record.START)
        self.assertEqual(25398279, record.END)
        self.assertEqual("KRAS", record.GENE)
        self.assertEqual("c.G40", record.CDS_MUTATION_SYNTAX)
        self.assertEqual("p.V14", record.AA_MUTATION_SYNTAX)
        self.assertEqual(ReportClass.region, record.REPORT)
        self.assertEqual("-", record.COMMENT)
        self.assertEqual("exon2", record.EXON)
        self.assertEqual("NM_004985", record.ACCESSION_NUMBER)
        self.assertEqual(1, len(record.TOTAL_DEPTH))
        self.assertEqual(["-"], record.TOTAL_DEPTH)

        record = reader.next()
        self.assertEqual("NC_000017.10", record.CHROMOSOME)
        self.assertEqual(37880986, record.START)
        self.assertEqual(37880987, record.END)
        self.assertEqual("ERBB2", record.GENE)
        self.assertEqual("-", record.CDS_MUTATION_SYNTAX)
        self.assertEqual("p.Y772", record.AA_MUTATION_SYNTAX)
        self.assertEqual(ReportClass.region_all, record.REPORT)
        self.assertEqual("-", record.COMMENT)
        self.assertEqual("exon20", record.EXON)
        self.assertEqual("NM_004448", record.ACCESSION_NUMBER)
        self.assertEqual(2, len(record.TOTAL_DEPTH))
        self.assertEqual(["-", "-"], record.TOTAL_DEPTH)

        with self.assertRaises(StopIteration):
            reader.next()
            55249003
#/data/ref_data/wp1/refFiles_20190130/refFiles/Mutations_Lung_20170125.csv:NC_000007.13	55242482	55242482	EGFR	c.C2252	p.T751	hotspot	-	exon19	NM_005228
#/data/ref_data/wp1/refFiles_20190130/refFiles/Mutations_Lung_20170125.csv:NC_000007.13	55242484	55242484	EGFR	c.T2254	p.S752	hotspot	-	exon19	NM_005228
#/data/ref_data/wp1/refFiles_20190130/refFiles/Mutations_Lung_20170125.csv:NC_000007.13	55242487	55242487	EGFR	c.C2257	p.P753	hotspot	-	exon19	NM_005228
#18-793	EGFR	frameshift substitution	exon20	-	c.2301_2336TGGTGCGGTCT	NM_005228	-	2-indel	yes	ok	775	763	12	0.0154838709677	-	-	-	No	-	-	-	-	-	-	-	-	-	-	-	-	+6|-6	-	-	NC_000007.13	55249003	5524903CAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGG	TGGTGCGGTCT	EGFR:NM_005228:exon20:c.2301_2336TGGTGCGGTCT

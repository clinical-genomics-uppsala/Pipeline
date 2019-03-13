# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import logging
import os
import shutil
import tempfile
import unittest

class TestPileup(unittest.TestCase):
    def setUp(self):
        # create fixtures
        self.tempdir = tempfile.mkdtemp()
        with open(os.path.join(self.tempdir,"multibp"),'w') as multibp:
            multibp.write("NC_000007.13\t55242465\t55242479\tGGAATTAAGAGAAGC\t-\tEGFR\tc.2235_2249del\tp.E746_A750del\tNM_005228.3\n")


    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_multibp_parsing(self):
         from src.lib.data.files.multibp import Reader
         reader = Reader(os.path.join(self.tempdir,"multibp"))
         multibp_variant = reader.next()

         self.assertEqual(multibp_variant.CHROMOSOME,"NC_000007.13")
         self.assertEqual(multibp_variant.START, "55242465")
         self.assertEqual(multibp_variant.STOP, "55242479")
         self.assertEqual(multibp_variant.REFERENCE, "GGAATTAAGAGAAGC")
         self.assertEqual(multibp_variant.VARIANT, "-")
         self.assertEqual(multibp_variant.GENE, "EGFR")
         self.assertEqual(multibp_variant.CDS_CHANGE, "c.2235_2249del")
         self.assertEqual(multibp_variant.AA_CHANGE, "p.E746_A750del")
         self.assertEqual(multibp_variant.TRANSCRIPT, "NM_005228")

    def test_multibp_data_structure(self):
          from src.lib.data.files.models import MultiBpVariantData
          data = MultiBpVariantData(os.path.join(self.tempdir,"multibp"))
          entry = data.get_data("NC_000007.13", "55242465", "55242479", "GGAATTAAGAGAAGC", "-")

          self.assertTrue(entry is not None)
          self.assertEqual(entry.CHROMOSOME,"NC_000007.13")
          self.assertEqual(entry.START, "55242465")
          self.assertEqual(entry.STOP, "55242479")
          self.assertEqual(entry.REFERENCE, "GGAATTAAGAGAAGC")
          self.assertEqual(entry.VARIANT, "-")
          self.assertEqual(entry.GENE, "EGFR")
          self.assertEqual(entry.CDS_CHANGE, "c.2235_2249del")
          self.assertEqual(entry.AA_CHANGE, "p.E746_A750del")
          self.assertEqual(entry.TRANSCRIPT, "NM_005228")

          entry = data.get_data("NC_000008.13", "55242465", "55242479", "GGAATTAAGAGAAGC", "-")
          self.assertTrue(entry is None)


if __name__ == '__main__':
    import logging
    import sys
    logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    unittest.main()

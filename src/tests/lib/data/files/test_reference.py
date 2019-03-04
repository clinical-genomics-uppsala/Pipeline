# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import logging
import os
import shutil
import tempfile
import unittest

logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

class TestReference(unittest.TestCase):
    def setUp(self):
        # create fixtures
        self.tempdir = tempfile.mkdtemp()
        with open(os.path.join(self.tempdir,"reference_info"),'w') as hotspots:
            hotspots.write("#Chr name\tNC\tID\tLength\n")
            hotspots.write("chr1\tNC_000001.9\tChr1#NC_000001.9#1#247249719#-1\t247249719\n")
            hotspots.write("chr2\tNC_000002.10\tChr2#NC_000002.10#1#242951149#-1\t242951149\n")
            hotspots.write("chr3\tNC_000003.10\tChr3#NC_000003.10#1#199501827#-1\t199501827\n")
            hotspots.write("chr4\tNC_000004.10\tChr4#NC_000004.10#1#191273063#-1\t191273063\n")
            hotspots.write("chr5\tNC_000005.8\tChr5#NC_000005.8#1#180857866#-1\t180857866\n")
            hotspots.write("chr6\tNC_000006.10\tChr6#NC_000006.10#1#170899992#-1\t170899992\n")
            hotspots.write("chr7\tNC_000007.12\tChr7#NC_000007.12#1#158821424#-1\t158821424")

    def tearDown(self):
        # delete fixtures
        shutil.rmtree(self.tempdir)

    def test_reference_info_import(self):
        from src.lib.data.files.reference import InfoImporter
        info = InfoImporter(os.path.join(self.tempdir,"reference_info"))
        self.assertTrue("chr1" in info)
        self.assertTrue("chr2" in info)
        self.assertTrue("chr3" in info)
        self.assertTrue("chr4" in info)
        self.assertTrue("chr5" in info)
        self.assertTrue("chr6" in info)
        self.assertTrue("chr7" in info)

        self.assertEqual(info['chr1']['length'], 247249719)
        self.assertEqual(info['chr4']['length'], 191273063)
        self.assertEqual(info['chr7']['length'], 158821424)

if __name__ == '__main__':
    import logging
    import sys
    logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    unittest.main()

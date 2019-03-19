# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import logging
import os
import shutil
import tempfile
import unittest
import gzip
class TestPileup(unittest.TestCase):
    def setUp(self):
        # create fixtures
        self.tempdir = tempfile.mkdtemp()
        with gzip.open(os.path.join(self.tempdir,"pileup.gz"),'wb') as hotspots:
            hotspots.write(b"chr1\t3628763\tG\t1293\t.......................................T....................................................................................................................................................................................................................,..................................,............................,.................,...............,,...........................................,,,,...........,,,,,.........,,,......,,,,,,,,,..,,,,,........,,,.........,.............,,,..........,,,,,,,,,,,,,......,,,,,...........,,,,,...............,,,,,,,..........,,,,,..............,,,,,,,...........,,.........,,,,,,,,,,,,................,,,,,,,,,,,,,,,..........,,.....................,,............,,,,,,,,,,.......,,,,........,,,,,,,,,,,,,...............,,,,,,,,,.......,,,,.....,,,,,,,,,,,,...........,,,,,,,,,,,,..........,,,,,,,,,,,,,,,...........,,,,,,,,,,,,,,........,,,,,,,,..,,,,..........,,,,,,,,,,........,,,,,,,.......,,,,,,,,,,,,,,.......................,,,,,,,,....,,,,,,,,,,,,,.....,,,,,,................,,,,,,,,,,,,,,,,,,,,,,,,.........................,,,,,,,,,,..........,,,,,,,,,,..................,,,,,,,,,,,...........,,,,,,,,,,,,,,,,,,...............,,,,,,,,,,,,,,.........,,,,,,,........................,,,,,,,,,,,,,,,,,,,....,,,,,,,,^].^].^].^].^].^].^].^].^].^].^].^].^].^],^],^],^],^],^],^],^],^],^],^],^],^],^],^],^],^],\t{NeN]z{{eNCCCe{eMNNCMN{{{NCCNNN@{e{]{{e!{N{{{{NNNNN@{e{{ee{{NNCMz{NNN{{{{eNP{{{NNNN{{NC@NNN{{ee{{bNNCNNNC{e{N{ONCCCP{{{{eeNNNCN{e{C{{N.Cz]{P{{{eee{NCy{NNCN{zNMCNCe{{{{PbebNNNNCCNCN;CCNN{NCNCNMNNCNNN.y{{{]NNCNCNe{{{{{NCCNCNNN{NNNNCCNNNNNNNNe{ee{e{eCNN{@!U{e{zbCNNN.N{eNCee{{{{NNNCNNNNCNCN!CzNCNC{eNNCNNNCCzbeNNNCNN.NCC{{NNCNCCNC{{e{NNN!{e{{zN{ee{{{NMN9!{{e{{{{{eNCCNNCCNNCNNNNNCNNe{N{{{{e{NCCNN.NCC!!b{e{NCNNNN.CNC!!{{zeNCCNCCNC{{{{bNCC!!!!!!!bNCN!!!{eeNNCCC@N!{{{{{z{NC@B{eNCNNNNMNNNC!!{eNMNNN.CNM.NNNC!!!!!!!y{{{NCC!!!!.{{{zeeNCCCCN!!!b{{eNNNNNCCNNNNN!!!!!!eeNCNCNNCCCCNC!{{{{eNCNMCNNCNNNNNN!!{NMCNN5CCNC!!z{{{{eeNNCNCCNNNC!!!!z{eNNNCCNNNNNNNCNN!!!!!!!!!!!!!{MNNNNNCNNNNe{ee{{ee{NNNNNCNNCNNN!!{eNNCCNNNNNCM!!!!!!!!!NCNCNCN!!!!{NNNCNNCN!!!!!!!!!!!!{{NNNNNCNC.NCNNNN!!!!!!!{eNNCNN!!!!{{eNNN!!!!!!!!!!!{{z{{{NNNNN@!!!!!!!!!!!{ee{NMNNMNNNN!@!!!!!!!!!!e{e{eNNNCCC!!!!!!!!!!!!!!z{z{{NNNN!!!!!!!eN!!!!NNNCCNNCCMNNNN!!!!!!{{eeNMNC@!!!!!!{CCNNCCNN.CC!!!!!!!!!e{{z{{yeWN.CNNNCCNNCNNCNNN!!!!!y{NN!!!!!!!!!!!!!ezNNN!!!!!!@Cze{e{{MCNNNMNCNN!!!!!!!!!!!!!!!!!!!!!!e{{{e{{CNNNNCNNNNNCCCCNNCC!!!!!!!!!{{{NNCNNNMNC!!!!!!!!{{e{eNNNCCCNNNCNCNMN!!!!!!!!!e{NNNCNNNN@NNN!!!!!!!!!!!!!!!z{bNN@NMM@NNN@N!!!!!!!!!!!!!!yzNNM@NNNN!!!!!!{{y_bNNMNN@NNN@NMM@NNNNN!!!!!!!!!!!!!!!!!!!b@NN!!!!!!!!{zb@N@N@NMNNNNCN!!!!!!!!!!!!!!\n")
            hotspots.write(b"chr1\t3628764\tA\t1316\t.$.$..........................................................................................................................................................................................................................................................,..........C.......................,............................,.................,...............,,...........................................,,,,...........,,,,,.........,,,......,,,,,,,,,..,,,,,........,,,.........,.............,,,..........,,,,,,,,,,,,,......,,,,,...........,,,c,...............,,,,,,,..........,,,,,..............,,,,,,,...........,,.........,,,,,,,,,,,,................,,,,,,,,,,,,,,,..........,,.....................,,............,,,,,,,,,,.......,,,,........,,,,,,,,,,,,,...............,,,,,,,,,.......,,,,.....,,,,,,,,,,,,...........,,,,,,,,,,,,..........,,,,t,,,,,,,,,,...........,,,,,,,,,,,,,,........,,,,,,,,..,,,,..........,,,,,,,,,,........,,,,,,,.......,,,,,,,,,,,,,,.......................,,,,,,,,....,,,,,,,,,,,,,.....,,,,,,................,,,,,,,,,,,,,,,,,,,,,,,,.........................,,,,,,,,,,..........,,,,,,,,,,..................,,,,,,,,,,,...........,,,,,,,,,,,,,,,,,,...............,,,,,,,,,,,,,,.........,,,,,,,........................,,,,,,,,,,,,,,,,,,,....,,,,,,,,.............,,,,,,,,,,,,,,,,,^].^].^].^].^].^].^].^].^].^].^],^],^],^],^],^],^],^],^],^],^],^],^],\t{{eNe{{{eNCCCe{eNNNCNN{{{NCCNNNC{e{;{{eP{N{{{{NNNNNC{e{{ee{{NNCW{zNNN{{{{eNe{{{{{{N{{NCCNNN{{]e{{eNN@NNNC{e{N{bNCCCe{{{{eeNNNCN{e{C{{NCC{e{<{{{eeb{NCz{NN@N{zNNCNCe{{{{_]ebNNNNCCNCNCCCNN{NCNCNMNNCNNN;z{{{b{{eNCNe{{{{{M@CNCNNN{NNNNCCNNNNNNNNe{eb{e{eCNN{5!e{e{{eCNNN.N{eNCee{{{{NNNCNNNNCNCN!CzNCN@{PNNCNNNCC{eeNNNCNNCNCC{{NNCNC@NC{{e{NNN!{e{{{Nzeb{{{NNN!!{{e{{{{{eNCCNNCCNNCNNNNNCNNe{N{{{{e{NCCNNCNCC!!e{e{NCNNNNCCNC!!{{{eNCCNCCNC{{{{;NCC!!!!!!!eNCN!!!{MeNNCCCCM!{{{{{z{NCC]{eN@NNNNNNNNC!!zeNNNNNCCNMCNNNC!!!!!!!z{{{NCC!!!!.{{{{ePNCCCCN!!!e{{eNNNNNCCNNNNN!!!!!!eeNCNCNNCCC@NC!{{{zeNCNMCNNCNNNNNN!!{NNCNNCCCNC!!{z{{{ebNNCNCCNNN@!!!!{{eNNNCCNNNNNNNCNN!!!!!!!!!!!!!{NNNNNNCNNNNe{ee{{ee{{NNNNCNNCNNN!!{eNNCCNNNNNCN!!!!!!!!!NCNCNCN!!!!{NNNCNNCN!!!!!!!!!!!!{{NNNNNCNC@NCMNNN!!!!!!!zeNN@NN!!!!{{eNNN!!!!!!!!!!!{{{{{{NNNNN;!!!!!!!!!!!{e]{NMNNNNNNM!.!!!!!!!!!!e{]{eNNNC@C!!!!!!!!!!!!!!z{z{{{NNN!!!!!!!eN!!!!NNNCCNNCCNNNNN!!!!!!{{eeNMNC@!!!!!!{CCNNCCNN..C!!!!!!!!!e{{{{{zePN.CNNNCCNNCNNCNNN!!!!!{{NN!!!!!!!!!!!!!eyNNN!!!!!!CCze{e{{NCNNNNNCNN!!!!!!!!!!!!!!!!!!!!!!e{{{e{{CNNNNCNNNNNCCCCNNC.!!!!!!!!!{{{NNCNNNMNC!!!!!!!!{{b{bNNNCCCNNNCNCNNN!!!!!!!!!e{NNNCNNNNCNNN!!!!!!!!!!!!!!!{{eNNCNNNCNNNCN!!!!!!!!!!!!!!yyNNM@NNNN!!!!!!{{zbbNNMNN@NNN@NMM@NNNNN!!!!!!!!!!!!!!!!!!!b@NN!!!!!!!!{z_bN5N@NMNNNNCN!!!!!!!!!!!!!!5yz{@NNNN@!!!!!!!!!!!!!\n")
            with open(os.path.join(self.tempdir,"pileup"),'w') as hotspots:
                hotspots.write("chr1\t3628763\tG\t1293\t.......................................T....................................................................................................................................................................................................................,..................................,............................,.................,...............,,...........................................,,,,...........,,,,,.........,,,......,,,,,,,,,..,,,,,........,,,.........,.............,,,..........,,,,,,,,,,,,,......,,,,,...........,,,,,...............,,,,,,,..........,,,,,..............,,,,,,,...........,,.........,,,,,,,,,,,,................,,,,,,,,,,,,,,,..........,,.....................,,............,,,,,,,,,,.......,,,,........,,,,,,,,,,,,,...............,,,,,,,,,.......,,,,.....,,,,,,,,,,,,...........,,,,,,,,,,,,..........,,,,,,,,,,,,,,,...........,,,,,,,,,,,,,,........,,,,,,,,..,,,,..........,,,,,,,,,,........,,,,,,,.......,,,,,,,,,,,,,,.......................,,,,,,,,....,,,,,,,,,,,,,.....,,,,,,................,,,,,,,,,,,,,,,,,,,,,,,,.........................,,,,,,,,,,..........,,,,,,,,,,..................,,,,,,,,,,,...........,,,,,,,,,,,,,,,,,,...............,,,,,,,,,,,,,,.........,,,,,,,........................,,,,,,,,,,,,,,,,,,,....,,,,,,,,^].^].^].^].^].^].^].^].^].^].^].^].^].^],^],^],^],^],^],^],^],^],^],^],^],^],^],^],^],^],\t{NeN]z{{eNCCCe{eMNNCMN{{{NCCNNN@{e{]{{e!{N{{{{NNNNN@{e{{ee{{NNCMz{NNN{{{{eNP{{{NNNN{{NC@NNN{{ee{{bNNCNNNC{e{N{ONCCCP{{{{eeNNNCN{e{C{{N.Cz]{P{{{eee{NCy{NNCN{zNMCNCe{{{{PbebNNNNCCNCN;CCNN{NCNCNMNNCNNN.y{{{]NNCNCNe{{{{{NCCNCNNN{NNNNCCNNNNNNNNe{ee{e{eCNN{@!U{e{zbCNNN.N{eNCee{{{{NNNCNNNNCNCN!CzNCNC{eNNCNNNCCzbeNNNCNN.NCC{{NNCNCCNC{{e{NNN!{e{{zN{ee{{{NMN9!{{e{{{{{eNCCNNCCNNCNNNNNCNNe{N{{{{e{NCCNN.NCC!!b{e{NCNNNN.CNC!!{{zeNCCNCCNC{{{{bNCC!!!!!!!bNCN!!!{eeNNCCC@N!{{{{{z{NC@B{eNCNNNNMNNNC!!{eNMNNN.CNM.NNNC!!!!!!!y{{{NCC!!!!.{{{zeeNCCCCN!!!b{{eNNNNNCCNNNNN!!!!!!eeNCNCNNCCCCNC!{{{{eNCNMCNNCNNNNNN!!{NMCNN5CCNC!!z{{{{eeNNCNCCNNNC!!!!z{eNNNCCNNNNNNNCNN!!!!!!!!!!!!!{MNNNNNCNNNNe{ee{{ee{NNNNNCNNCNNN!!{eNNCCNNNNNCM!!!!!!!!!NCNCNCN!!!!{NNNCNNCN!!!!!!!!!!!!{{NNNNNCNC.NCNNNN!!!!!!!{eNNCNN!!!!{{eNNN!!!!!!!!!!!{{z{{{NNNNN@!!!!!!!!!!!{ee{NMNNMNNNN!@!!!!!!!!!!e{e{eNNNCCC!!!!!!!!!!!!!!z{z{{NNNN!!!!!!!eN!!!!NNNCCNNCCMNNNN!!!!!!{{eeNMNC@!!!!!!{CCNNCCNN.CC!!!!!!!!!e{{z{{yeWN.CNNNCCNNCNNCNNN!!!!!y{NN!!!!!!!!!!!!!ezNNN!!!!!!@Cze{e{{MCNNNMNCNN!!!!!!!!!!!!!!!!!!!!!!e{{{e{{CNNNNCNNNNNCCCCNNCC!!!!!!!!!{{{NNCNNNMNC!!!!!!!!{{e{eNNNCCCNNNCNCNMN!!!!!!!!!e{NNNCNNNN@NNN!!!!!!!!!!!!!!!z{bNN@NMM@NNN@N!!!!!!!!!!!!!!yzNNM@NNNN!!!!!!{{y_bNNMNN@NNN@NMM@NNNNN!!!!!!!!!!!!!!!!!!!b@NN!!!!!!!!{zb@N@N@NMNNNNCN!!!!!!!!!!!!!!\n")
                hotspots.write("chr1\t3628764\tA\t1316\t.$.$..........................................................................................................................................................................................................................................................,..........C.......................,............................,.................,...............,,...........................................,,,,...........,,,,,.........,,,......,,,,,,,,,..,,,,,........,,,.........,.............,,,..........,,,,,,,,,,,,,......,,,,,...........,,,c,...............,,,,,,,..........,,,,,..............,,,,,,,...........,,.........,,,,,,,,,,,,................,,,,,,,,,,,,,,,..........,,.....................,,............,,,,,,,,,,.......,,,,........,,,,,,,,,,,,,...............,,,,,,,,,.......,,,,.....,,,,,,,,,,,,...........,,,,,,,,,,,,..........,,,,t,,,,,,,,,,...........,,,,,,,,,,,,,,........,,,,,,,,..,,,,..........,,,,,,,,,,........,,,,,,,.......,,,,,,,,,,,,,,.......................,,,,,,,,....,,,,,,,,,,,,,.....,,,,,,................,,,,,,,,,,,,,,,,,,,,,,,,.........................,,,,,,,,,,..........,,,,,,,,,,..................,,,,,,,,,,,...........,,,,,,,,,,,,,,,,,,...............,,,,,,,,,,,,,,.........,,,,,,,........................,,,,,,,,,,,,,,,,,,,....,,,,,,,,.............,,,,,,,,,,,,,,,,,^].^].^].^].^].^].^].^].^].^].^],^],^],^],^],^],^],^],^],^],^],^],^],\t{{eNe{{{eNCCCe{eNNNCNN{{{NCCNNNC{e{;{{eP{N{{{{NNNNNC{e{{ee{{NNCW{zNNN{{{{eNe{{{{{{N{{NCCNNN{{]e{{eNN@NNNC{e{N{bNCCCe{{{{eeNNNCN{e{C{{NCC{e{<{{{eeb{NCz{NN@N{zNNCNCe{{{{_]ebNNNNCCNCNCCCNN{NCNCNMNNCNNN;z{{{b{{eNCNe{{{{{M@CNCNNN{NNNNCCNNNNNNNNe{eb{e{eCNN{5!e{e{{eCNNN.N{eNCee{{{{NNNCNNNNCNCN!CzNCN@{PNNCNNNCC{eeNNNCNNCNCC{{NNCNC@NC{{e{NNN!{e{{{Nzeb{{{NNN!!{{e{{{{{eNCCNNCCNNCNNNNNCNNe{N{{{{e{NCCNNCNCC!!e{e{NCNNNNCCNC!!{{{eNCCNCCNC{{{{;NCC!!!!!!!eNCN!!!{MeNNCCCCM!{{{{{z{NCC]{eN@NNNNNNNNC!!zeNNNNNCCNMCNNNC!!!!!!!z{{{NCC!!!!.{{{{ePNCCCCN!!!e{{eNNNNNCCNNNNN!!!!!!eeNCNCNNCCC@NC!{{{zeNCNMCNNCNNNNNN!!{NNCNNCCCNC!!{z{{{ebNNCNCCNNN@!!!!{{eNNNCCNNNNNNNCNN!!!!!!!!!!!!!{NNNNNNCNNNNe{ee{{ee{{NNNNCNNCNNN!!{eNNCCNNNNNCN!!!!!!!!!NCNCNCN!!!!{NNNCNNCN!!!!!!!!!!!!{{NNNNNCNC@NCMNNN!!!!!!!zeNN@NN!!!!{{eNNN!!!!!!!!!!!{{{{{{NNNNN;!!!!!!!!!!!{e]{NMNNNNNNM!.!!!!!!!!!!e{]{eNNNC@C!!!!!!!!!!!!!!z{z{{{NNN!!!!!!!eN!!!!NNNCCNNCCNNNNN!!!!!!{{eeNMNC@!!!!!!{CCNNCCNN..C!!!!!!!!!e{{{{{zePN.CNNNCCNNCNNCNNN!!!!!{{NN!!!!!!!!!!!!!eyNNN!!!!!!CCze{e{{NCNNNNNCNN!!!!!!!!!!!!!!!!!!!!!!e{{{e{{CNNNNCNNNNNCCCCNNC.!!!!!!!!!{{{NNCNNNMNC!!!!!!!!{{b{bNNNCCCNNNCNCNNN!!!!!!!!!e{NNNCNNNNCNNN!!!!!!!!!!!!!!!{{eNNCNNNCNNNCN!!!!!!!!!!!!!!yyNNM@NNNN!!!!!!{{zbbNNMNN@NNN@NMM@NNNNN!!!!!!!!!!!!!!!!!!!b@NN!!!!!!!!{z_bN5N@NMNNNNCN!!!!!!!!!!!!!!5yz{@NNNN@!!!!!!!!!!!!!\n")

    def tearDown(self):
        shutil.rmtree(self.tempdir)


    def test_header_parsing(self):
         from src.lib.data.files.pileup import Reader
         reader = Reader(os.path.join(self.tempdir,"pileup"))
         pileup_position = reader.next()

         self.assertEqual(pileup_position.CHROMOSOME, "chr1")
         self.assertEqual(pileup_position.POSITION, 3628763)
         self.assertEqual(pileup_position.REFERENCE_BASE, "G")
         self.assertEqual(pileup_position.DEPTH, 1293)

         pileup_position = reader.next()
         self.assertEqual(pileup_position.CHROMOSOME, "chr1")
         self.assertEqual(pileup_position.POSITION, 3628764)
         self.assertEqual(pileup_position.REFERENCE_BASE, "A")
         self.assertEqual(pileup_position.DEPTH, 1316)

         with self.assertRaises(StopIteration):
             reader.next()

    def test_header_parsing_gz(self):
         from src.lib.data.files.pileup import Reader
         reader = Reader(os.path.join(self.tempdir,"pileup.gz"))
         pileup_position = reader.next()

         self.assertEqual(pileup_position.CHROMOSOME, "chr1")
         self.assertEqual(pileup_position.POSITION, 3628763)
         self.assertEqual(pileup_position.REFERENCE_BASE, "G")
         self.assertEqual(pileup_position.DEPTH, 1293)

         pileup_position = reader.next()
         self.assertEqual(pileup_position.CHROMOSOME, "chr1")
         self.assertEqual(pileup_position.POSITION, 3628764)
         self.assertEqual(pileup_position.REFERENCE_BASE, "A")
         self.assertEqual(pileup_position.DEPTH, 1316)

         with self.assertRaises(StopIteration):
             reader.next()


if __name__ == '__main__':
    import logging
    import sys
    logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    unittest.main()

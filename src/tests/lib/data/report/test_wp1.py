# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import logging
import os
import shutil
import tempfile
import unittest

from pysam import tabix_index
logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

class TestWp1Reports(unittest.TestCase):
    def setUp(self):
        # create fixtures
        self.tempdir = tempfile.mkdtemp()

        self.multibp = os.path.join(self.tempdir,"multibp")
        with open(os.path.join(self.tempdir,"multibp"),'w') as multibp:
            multibp.write("NC_000007.13\t55242465\t55242479\tGGAATTAAGAGAAGC\t-\tEGFR\tc.2235_2249del\tp.E746_A750del\tNM_005228.3\n")

        self.hotspot = os.path.join(self.tempdir,"hotspot")
        with open(self.hotspot,'w') as hotspots:
            hotspots.write("#Chr\tStart\tEnd\tGene\tCDS_mutation_syntax\tAA_mutation_syntax\tReport\tcomment\tExon\tAccession_number\n")
            hotspots.write("NC_000002.11\t29445258\t29445258\tALK\tc.G3467\tp.C1156\tregion\tresistance_mutation\texon22\tNM_004304\n")
            hotspots.write("NC_000007.13\t55242466\t55242466\tEGFR\tc.G2236\tp.E746\thotspot\t-\texon19\tNM_005228\n")
            hotspots.write("NC_000007.13\t140453136\t140453136\tBRAF\tc.T1799\tp.V600\thotspot\t-\texon15\tNM_004333\n")
            hotspots.write("NC_000007.13\t116412043\t116412043\tMET\tc.G3082\tp.D1028\thotspot\t-\texon14\tNM_001127500\n")

        self.pileup = os.path.join(self.tempdir,"pilup")
        with open(self.pileup, "w") as pileup:
            pileup.write("chr2\t29445216\tC\t601\t...\t...\n")
            pileup.write("chr2\t29445258\tC\t600\t...\t...\n")
            pileup.write("chr7\t55242464\tC\t830\t...\t...\n")
            pileup.write("chr7\t55242465\tC\t831\t...\t...\n")
            pileup.write("chr7\t55242466\tC\t832\t...\t...\n")
            pileup.write("chr7\t55242467\tC\t833\t...\t...\n")
            pileup.write("chr7\t55242468\tC\t834\t...\t...\n")
            pileup.write("chr7\t55242469\tC\t835\t...\t...\n")
            pileup.write("chr7\t55242470\tC\t836\t...\t...\n")
            pileup.write("chr7\t55242471\tC\t837\t...\t...\n")
            pileup.write("chr7\t55242472\tC\t838\t...\t...\n")
            pileup.write("chr7\t55242473\tC\t839\t...\t...\n")
            pileup.write("chr7\t140453136\tA\t834\t...\t...\n")
            pileup.write("chr7\t116412043\tA\t999\t...\t...\n")
            pileup.write("chr7\t116411979\tT\t998\t...\t...\n")
            pileup.write("chr17\t37883141\tG\t1288\t.$.$,$,$.....,,,...,,,,,,,,..,,,..................,,,,,,.............,,......,,,,,,,,....,,,,,,,,....,,,,,,,,,,,,,,............,,,,,,,,,.....,,,......,,,...,,,.......,,,,,,,,,,,......,,,,,,,,,,....,,,,...,,,,,,,,,,...,,,.$.$........,,,,,$,$,.......,,,,,,,,.............,,,,,,,,,.,.......,,,,t,,......,,,,,,,,..................,,,,,,,,,,,,,,,,,,...,,.................,,,,,,,,.,,,,,,,,.........,,,,,,,,..........,,,,,,,,,,,,,.....,,,,,,,,,.............,...,,.......,,,,,,..........,,,,,..,.,,,,,,,,,,........,,,,,,,,.....,,,,,,,,,,.,.........,,,,,,,,,,............,,,,,,,,................,,,,,,,,,,......,,,,,,........,,,,,,,.......,,,,,,,,,,......,,,,,,,.........,,,,,,,,,,,,...........,,,,,,,,,,,,,,,,,,.....,,,,,,,,........,,,,,,,,...,,,,,,,,,,,,........,,,,,,,,,,,..,,,,,,,..,,,,...........,,,,,........,,,,...,,,,.......,,,,,,,,..........,,,,,,,,,......,,,,,,,........,,,,.........,,,,,,,,,.....,,,,,........,,,,,,,,......,,,,,,,,,.....,,,,,,,........,,,,,,,,,,....,,,..........,,,,,,,,...........,,,,,,,,,,,,,,,,.....,,,,,,.......,,.....,,,,,,,....,........,,,,,,,,,......,,,,,,,,,.....,,,,,,,,,.........,,,,,,,,,,.....,,,....,,,,....,,,,,,,.....,,,,,,........,,,,,,,,,,,,.....,,......,,,,,,,,,...,,,,,,,,,,,,,......,,,,..,,,,,..,,,,,,,.,,,,,,,.,,,,,,,,,,....,,,,,.,,,.,,,,,...,,,,,,,,^].^],^],^],^],^],^],\t{eNN{e{e{NNN{{eNNNNMNNM{NNNNbz{{ee{ze{{zee{{{eNNNN@!{{e{{{{{CNNNNCN{{{e{{NNNNNNN@{eezNNCCCNNNb{e{NNCNNNCCNNNNC!{{{{{{{zz{NNMNN!!!!!C{{e{{CNC{{{{CNN!!{eeNN!{{{{{b{N@CNN@NCNN!{eeeb{NNNCNNNN!!{{{eNNN!{{{NCNNN!!!!!e{NNNCzzee{{{eeNNCCM!!!{e{{bNCNNNCNN!!{{e{{{{NCCNNCNCCC!!!!!{Ny{{{{NNNCNN.C!ee{{{NNNMMN!!!{e{{Pe{{{{z{{z{zNNNNCNNCCN!!!!!!!!!!{{yN!e{{ee{{{{e{{{{eNMNNNCN!!!eNNN@NN!!{{{{{N@CCN!!!!!!!z{{{{eMNCNNNNNNCMNC!!!!{{{{eC!!!!!!!!{{e{{{{NMCCCNNb{{N!{{e{{eCN!!!!!{{{e{NNNCCN!!!!NC!bC!!!!!!!!!{{{{{NNNNN!!!!!!e{e{eN!!!!!!!!!{!{{e{{{{bNNCNC!!!!!!b{e{e{{{{{xe!!!!!!!!{b{{{{e{e{{{{eNNNNNNN!!!!!e{{{{eN!!!!!zee{{eNNCN!!!!!{{{{NCCNNNC!!!!!!{{zeeeNCNN!!!e{{{{e{{CNNNN!!!!!!!!{ze{{{{{@NNNNCNNCCNCN!!!!!!!!{e{NCNNN!!!!!eeNNNCNCN!!!!!!!{e{N!!!!!!!!!!!ee{{{bNNNNNNNC!!!!!e{N!!!!!!{N!!!!{{{{z{{{eNN!!!!!{{e{{{NN!!!!{{eN!!!ee{{NNNNNCN!!!!e{{NNNNNNCNNCNN!!!!{{{{NNN!!!!!!{{{NNCCNN!!!{{{NCNCNMNN!!!!!!!e{e{CN!!!!{be{{NNNN!!!!!!!{{e{NCN.C!!!!!!ezNCC!!!!!!!e{{{{{zbN@!!!!!!!!{NCN!!!P{{e{NNNNNNNCC!!!!eee{{zNNNNNNCCNN!!!!!!!!!!!{{NNNC!!!!!eNCNNNC!!e{eeNNN!!!!!{NNC!b{{eNCNCNN!!!!!!!yeNNCM!!!!!!!!!{NNNN!!!!!N!!!z{eezeNNNNN!!!!!!!!{e{NN!!!eeee!!!!z{NNN!!!!!!eeNNC!!!!!!e{{eeNNN!!!!!!!!!!!!ezNCC!!ee{ePC!!!!!!!!!eCN!!!!!!!!!!!!!{eNCNC!!!!ez!!!!!NNNNNCNN!C!!!!!!!CNCC@NM@!!!bbM@!!!!!@!!!Z!!!!!@@@NCN!!!!!bMN!!!!\n")
        self.vcf = os.path.join(self.tempdir,"vcfdata.vcf")
        with open(self.vcf, "w", encoding="ascii") as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write("##reference=hg19\n")
            vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcf.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n")
            vcf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n")
            vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n")
            vcf.write("##contig=<ID=chr1,length=249250621,assembly=hg19>\n")
            vcf.write("##contig=<ID=chr2,length=242951149,assembly=hg19>\n")
            vcf.write("##contig=<ID=chr3,length=199501827,assembly=hg19>\n")
            vcf.write("##contig=<ID=chr7,length=158821424,assembly=hg19>\n")
            vcf.write("##contig=<ID=chr12,length=133851895,assembly=hg19>\n")
            vcf.write("##contig=<ID=chr17,length=81195210,assembly=hg19>\n")
            vcf.write("##INFO=<ID=ANNOVAR_DATE,Number=1,Type=String,Description=\"Flag the start of ANNOVAR annotation for one alternative allele\">\n")
            vcf.write("##INFO=<ID=Func.refGene,Number=.,Type=String,Description=\"Func.refGene annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=Gene.refGene,Number=.,Type=String,Description=\"Gene.refGene annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=GeneDetail.refGene,Number=.,Type=String,Description=\"GeneDetail.refGene annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=ExonicFunc.refGene,Number=.,Type=String,Description=\"ExonicFunc.refGene annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=AAChange.refGene,Number=.,Type=String,Description=\"AAChange.refGene annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=1000g2015aug_eur,Number=1,Type=Float,Description=\"1000g2015aug_eur annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=snp138,Number=.,Type=String,Description=\"snp138 annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=snp138NonFlagged,Number=.,Type=String,Description=\"snp138NonFlagged annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=esp6500siv2_ea,Number=1,Type=Float,Description=\"esp6500siv2_ea annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=cosmic70,Number=.,Type=String,Description=\"cosmic70 annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=clinvar_20150629,Number=.,Type=String,Description=\"clinvar_20150629 annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=ALLELE_END,Number=0,Type=Flag,Description=\"Flag the end of ANNOVAR annotation for one alternative allele\">\n")
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
            vcf.write(r"chr2	29445216	.	ATCAGGGCTTCCATGAGGAAATCCAGTTCGTCCTGTTCAGAGCACACTTCAGGCAGCGTCTGGGCAGAGAAGGGGAGGGTGGGGAGGAGGAG	GCTGGTGCGGTCTCA	50	.	ANNOVAR_DATE=2015-06-17;Func.refGene=exonic;Gene.refGene=ALK;GeneDetail.refGene=.;ExonicFunc.refGene=frameshift_substitution;AAChange.refGene=ALK:NM_004304:exon22:c.3451_3509TGAGACCGCACCAGC;1000g2015aug_eur=.;snp138=.;snp138NonFlagged=.;esp6500siv2_ea=.;cosmic70=.;clinvar_20150629=.;ALLELE_END	GT:DP:AD	0/1:596:580,16")
            vcf.write("\n")
            vcf.write(r"chr7	55242464	rs121913421;COSM6223	AGGAATTAAGAGAAGC	A	60	.	ANNOVAR_DATE=2015-06-17;Func.refGene=exonic;Gene.refGene=EGFR;GeneDetail.refGene=.;ExonicFunc.refGene=nonframeshift_deletion;AAChange.refGene=EGFR:NM_005228:exon19:c.2235_2249del:p.745_750del;1000g2015aug_eur=.;snp138=rs121913421;snp138NonFlagged=.;esp6500siv2_ea=.;cosmic70=ID\x3dCOSM6223\x3bOCCURENCE\x3d894(lung),2(thyroid),3(upper_aerodigestive_tract),1(oesophagus),2(ovary),4(salivary_gland);clinvar_20150629=CLINSIG\x3ddrug-response\x3bCLNDBN\x3dTyrosine_kinase_inhibitor_response\x3bCLNREVSTAT\x3dno_assertion_criteria_provided\x3bCLNACC\x3dRCV000150617.1\x3bCLNDSDB\x3dMedGen\x3bCLNDSDBID\x3dCN225347;ALLELE_END	GT:DP:AD	0/1:613:214,399")
            vcf.write("\n")
            vcf.write(r"chr7	116411979	.	TGTAAGCCCAACTACAGAAATGGTTTCAAATGAATCTGTAGACTACCGAGCTACTTTTCCAGAAGGTATATTTCAGTTTATTGTTCTGAGAA	CT	70	.	ANNOVAR_DATE=2015-06-17;Func.refGene=exonic;Gene.refGene=MET;GeneDetail.refGene=.;ExonicFunc.refGene=nonframeshift_substitution;AAChange.refGene=MET:NM_000245:exon14:c.2964_3028CT,MET:NM_001127500:exon14:c.3018_3082CT;1000g2015aug_eur=.;snp138=.;snp138NonFlagged=.;esp6500siv2_ea=.;cosmic70=.;clinvar_20150629=.;ALLELE_END	GT:DP:AD	0/1:391:386,5")
            vcf.write("\n")
            vcf.write(r"chr7	140453136	rs113488022;COSM476	A	T	80	.	.;ANNOVAR_DATE=2015-06-17;Func.refGene=exonic;Gene.refGene=BRAF;GeneDetail.refGene=.;ExonicFunc.refGene=nonsynonymous_SNV;AAChange.refGene=BRAF:NM_004333:exon15:c.T1799A:p.V600E;1000g2015aug_eur=.;snp138=rs113488022;snp138NonFlagged=.;esp6500siv2_ea=.;cosmic70=ID\x3dCOSM476\x3bOCCURENCE\x3d3(adrenal_gland),3(urinary_tract),4(genital_tract),7(pancreas),16(liver),68(eye),2(prostate),4161(large_intestine),21(biliary_tract),11(breast),4377(skin),9(small_intestine),1(autonomic_ganglia),20(pituitary),242(ovary),32(soft_tissue),7(testis),459(haematopoietic_and_lymphoid_tissue),9534(thyroid),78(lung),1(meninges),12(upper_aerodigestive_tract),717(NS),15(bone),2(oesophagus),4(stomach),228(central_nervous_system);clinvar_20150629=CLINSIG\x3dpathogenic,pathogenic|pathogenic|pathogenic|pathogenic|pathogenic|pathogenic|pathogenic\x3bCLNDBN\x3dRasopathy,Carcinoma_of_colon|Papillary_thyroid_carcinoma|Astrocytoma\x2c_low-grade\x2c_somatic|Germ_cell_tumor\x2c_nonseminomatous|Non-small_cell_lung_cancer|Malignant_melanoma|not_provided\x3bCLNREVSTAT\x3dcriteria_provided\x2c_single_submitter,no_assertion_criteria_provided|no_assertion_criteria_provided|no_assertion_criteria_provided|no_assertion_criteria_provided|no_assertion_criteria_provided|no_assertion_criteria_provided|criteria_provided\x2c_single_submitter\x3bCLNACC\x3dRCV000033335.2,RCV000014992.10|RCV000014993.10|RCV000014994.10|RCV000022677.10|RCV000037936.2|RCV000067669.8|RCV000080903.3\x3bCLNDSDB\x3dMedGen,MedGen:SNOMED_CT|MedGen:OMIM:Orphanet:SNOMED_CT|.|MedGen|MedGen:SNOMED_CT|MedGen:SNOMED_CT|MedGen\x3bCLNDSDBID\x3dCN166718,C0699790:269533000|C0238463:188550:ORPHA146:255029007|.|C1266158|C0007131:254637007|C0025202:2092003|CN221809;ALLELE_END	GT:DP:AD	0/1:6142:4397,1744")
            vcf.write("\n")
            vcf.write(r"chr17	37883141	rs140272156	G	A	90	.	.;ANNOVAR_DATE=2015-06-17;Func.refGene=exonic;Gene.refGene=ERBB2;GeneDetail.refGene=.;ExonicFunc.refGene=nonsynonymous_SNV;AAChange.refGene=ERBB2:NM_001289937:exon25:c.G3044A:p.G1015E,ERBB2:NM_004448:exon25:c.G3044A:p.G1015E,ERBB2:NM_001005862:exon28:c.G2954A:p.G985E,ERBB2:NM_001289936:exon29:c.G2999A:p.G1000E;1000g2015aug_eur=0.004;snp138=rs140272156;snp138NonFlagged=rs140272156;esp6500siv2_ea=0.0020;cosmic70=.;clinvar_20150629=.;ALLELE_END	GT:DP:AD	0/1:3797:1824,1971")

        tabix_index(self.vcf , preset="vcf")
        self.reference = os.path.join(self.tempdir,"reference_info")
        with open(self.reference,'w') as reference:
            reference.write("#Chr name\tNC\tID\tLength\n")
            reference.write("chr1\tNC_000001.10\tchr1#NC_000001.10#1#249250621#-1\t249250621\n")
            reference.write("chr2\tNC_000002.11\tChr2#NC_000002.10#1#242951149#-1\t242951149\n")
            reference.write("chr3\tNC_000003.10\tChr3#NC_000003.10#1#199501827#-1\t199501827\n")
            reference.write("chr4\tNC_000004.10\tChr4#NC_000004.10#1#191273063#-1\t191273063\n")
            reference.write("chr5\tNC_000005.8\tChr5#NC_000005.8#1#180857866#-1\t180857866\n")
            reference.write("chr6\tNC_000006.10\tChr6#NC_000006.10#1#170899992#-1\t170899992\n")
            reference.write("chr7\tNC_000007.13\tchr7#NC_000007.13#1#159138663#-1\t159138663\n")
            reference.write("chr12\tNC_000012.11\tchr12#NC_000012.11#1#133851895#-1\t133851895\n")
            reference.write("chr17\tNC_000017.10\tchr17#NC_000017.10#1#81195210#-1\t81195210");

        self.hotspot_2 = os.path.join(self.tempdir,"hotspot_2")
        with open(self.hotspot_2,'w') as hotspots:
            hotspots.write("#Chr\tStart\tEnd\tGene\tCDS_mutation_syntax\tAA_mutation_syntax\tReport\tcomment\tExon\tAccession_number\n")
            hotspots.write("NC_000012.11\t25380168\t25380346\tKRAS\t-\t-\tindel\t-\texon3\tNM_004985\n")

        self.vcf_2 = os.path.join(self.tempdir,"vcf_2")
        with open(self.vcf_2,'w') as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write("##reference=hg19\n")
            vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcf.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n")
            vcf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n")
            vcf.write("##contig=<ID=chr12,length=133851895,assembly=hg19>\n")
            vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n")
            vcf.write("##INFO=<ID=ANNOVAR_DATE,Number=1,Type=String,Description=\"Flag the start of ANNOVAR annotation for one alternative allele\">\n")
            vcf.write("##INFO=<ID=Func.refGene,Number=.,Type=String,Description=\"Func.refGene annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=Gene.refGene,Number=.,Type=String,Description=\"Gene.refGene annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=GeneDetail.refGene,Number=.,Type=String,Description=\"GeneDetail.refGene annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=ExonicFunc.refGene,Number=.,Type=String,Description=\"ExonicFunc.refGene annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=AAChange.refGene,Number=.,Type=String,Description=\"AAChange.refGene annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=1000g2015aug_eur,Number=1,Type=Float,Description=\"1000g2015aug_eur annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=snp138,Number=.,Type=String,Description=\"snp138 annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=snp138NonFlagged,Number=.,Type=String,Description=\"snp138NonFlagged annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=esp6500siv2_ea,Number=1,Type=Float,Description=\"esp6500siv2_ea annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=cosmic70,Number=.,Type=String,Description=\"cosmic70 annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=clinvar_20150629,Number=.,Type=String,Description=\"clinvar_20150629 annotation provided by ANNOVAR\">\n")
            vcf.write("##INFO=<ID=ALLELE_END,Number=0,Type=Flag,Description=\"Flag the end of ANNOVAR annotation for one alternative allele\">\n")
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
            vcf.write(r"chr12	25380316	.	C	A	54	PASS	DP=4032;AF=0.001984;SB=21;DP4=2016,2008,0,8;AD=4024,8;ANNOVAR_DATE=2015-06-17;Func.refGene=exonic;Gene.refGene=KRAS;GeneDetail.refGene=.;ExonicFunc.refGene=stopgain;AAChange.refGene=KRAS:NM_004985:exon3:c.G142T:p.G48X,KRAS:NM_033360:exon3:c.G142T:p.G48X;1000g2015aug_eur=.;snp138=.;snp138NonFlagged=.;esp6500siv2_ea=.;cosmic70=.;clinvar_20150629=.;ALLELE_END	GT:AF:DP4:DP:AD	0/1:0.001984:2016,2008,0,8:4032:4024,8")
        tabix_index(self.vcf_2 , preset="vcf")

        self.pileup_2 = os.path.join(self.tempdir,"pileup_2")
        with open(os.path.join(self.tempdir,"pileup_2"),'w') as pileup:
            pileup.write("chr12\t25380315\tC\t4031\t.\t.\n")
            pileup.write("chr12\t25380316\tC\t4032\t.\t.\n")
            pileup.write("chr12\t25380317\tC\t4033\t.\t.")

    def tearDown(self):
        # delete fixtures
        shutil.rmtree(self.tempdir)

    def test_get_read_level(self):
        from src.lib.data.report.wp1 import get_read_level
        levels = [(300 ,"ok","yes"), (30 ,"low","yes"),(0 ,"low","not analyzable")]

        self.assertEqual(get_read_level(levels, 300), ("ok","yes"))
        self.assertEqual(get_read_level(levels, 200), ("low","yes"))
        self.assertEqual(get_read_level(levels, 30), ("low","yes"))
        self.assertEqual(get_read_level(levels, 15), ("low","not analyzable"))
        self.assertEqual(get_read_level(levels, 0), ("low","not analyzable"))
        self.assertEqual(get_read_level(levels, -15), ("-","zero"))
        self.assertEqual(get_read_level(levels, "-"), ("-","zero"))


    def test_filtered_mutation_creation_annovar(self):
        from src.lib.data.report.wp1 import generate_filtered_mutations
        levels = [(300 ,"ok","yes"), (30 ,"low","yes"),(0 ,"low","not analyzable")]

        report = os.path.join(self.tempdir,"filtered.report")
        generate_filtered_mutations("sample1", "variant_caller", report, levels, self.hotspot, self.vcf + ".gz" , self.pileup, self.reference, self.multibp, {"ERBB2": "NM_004448", "MET": "NM_001127500"})

        self.maxDiff = 10000
        with open(report, 'r') as report_result:
            head = report_result.readline()
            self.assertEqual(head.rstrip(), "\t".join(["Sample","Gene","Variant_type","Exon","AA_change",
                                                    "CDS_change","Accession_number","Comment","Report",
                                                    "Found","Min_read_depth300","Pileup_depth","Total_read_depth",
                                                    "Reference_read_depth","Variant_read_depth","Variant_allele_ratio", "Qual", "Caller",
                                                    "dbSNP_id","Ratio_1000G","Ratio_ESP6500","Clinically_flagged_dbSNP",
                                                    "Cosmic","ClinVar_CLNDB","Clinval_CLINSIG","Reference_plus_amplicons",
                                                    "Reference_minus_amplicons","Variant_plus_amplicons","Variant_minus_amplicons",
                                                    "Strands_A_F+F-S+S-","Strands_G_F+F-S+S-","Strands_C_F+F-S+S-","Strands_T_F+F-S+S-",
                                                    "Strands_Ins","Strands_Del","Ref_aligned_amplicons","Var_aligned_amplicons",
                                                    "Chr","Start","End","Reference_base","Variant_base","All_transcripts_annotation"]))

            result = report_result.readlines()
            self.assertEqual(len(result),5)
            ## ToDo make sure that depth is mean of whole deletion?
            self.assertEqual(result[0].rstrip(),
                "sample1\tEGFR\tnonframeshift_deletion\texon19\tp.E746_A750del\tc.2235_2249del\tNM_005228\t-\t2-indel\tyes\tok\t830\t613\t214\t399\t0.6508972267536705\t60.0\tvariant_caller\trs121913421\t-\t-\tYes\tID=COSM6223;OCCURENCE=894(lung),2(thyroid),3(upper_aerodigestive_tract),1(oesophagus),2(ovary),4(salivary_gland)\tTyrosine_kinase_inhibitor_response\tdrug-response\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNC_000007.13\t55242465\t55242479\tGGAATTAAGAGAAGC\t-\tEGFR:NM_005228:exon19:c.2235_2249del:p.745_750del")
            self.assertEqual(result[1].rstrip(),
                "sample1\tBRAF\tnonsynonymous_SNV\texon15\tp.V600E\tc.T1799A\tNM_004333\t-\t1-hotspot\tyes\tok\t834\t6142\t4397\t1744\t0.28394659719960924\t80.0\tvariant_caller\trs113488022\t-\t-\tYes\tID=COSM476;OCCURENCE=3(adrenal_gland),3(urinary_tract),4(genital_tract),7(pancreas),16(liver),68(eye),2(prostate),4161(large_intestine),21(biliary_tract),11(breast),4377(skin),9(small_intestine),1(autonomic_ganglia),20(pituitary),242(ovary),32(soft_tissue),7(testis),459(haematopoietic_and_lymphoid_tissue),9534(thyroid),78(lung),1(meninges),12(upper_aerodigestive_tract),717(NS),15(bone),2(oesophagus),4(stomach),228(central_nervous_system)\tRasopathy,Carcinoma_of_colon|Papillary_thyroid_carcinoma|Astrocytoma\x2c_low-grade\x2c_somatic|Germ_cell_tumor\x2c_nonseminomatous|Non-small_cell_lung_cancer|Malignant_melanoma|not_provided\tpathogenic,pathogenic|pathogenic|pathogenic|pathogenic|pathogenic|pathogenic|pathogenic\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNC_000007.13\t140453136\t140453136\tA\tT\tBRAF:NM_004333:exon15:c.T1799A:p.V600E")
            self.assertEqual(result[2].rstrip(),
                "sample1\tMET\tnonframeshift_substitution\texon14\t-\tc.3018_3082CT\tNM_001127500\t-\t2-indel\tyes\tok\t998\t391\t386\t5\t0.01278772378516624\t70.0\tvariant_caller\t-\t-\t-\tNo\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNC_000007.13\t116411979\t116412070\tTGTAAGCCCAACTACAGAAATGGTTTCAAATGAATCTGTAGACTACCGAGCTACTTTTCCAGAAGGTATATTTCAGTTTATTGTTCTGAGAA\tCT\tMET:NM_000245:exon14:c.2964_3028CT,MET:NM_001127500:exon14:c.3018_3082CT")
            self.assertEqual(result[3].rstrip(),
                "sample1\tALK\tframeshift_substitution\texon22\t-\tc.3451_3509TGAGACCGCACCAGC\tNM_004304\tresistance_mutation\t3-check\tyes\tok\t601\t596\t580\t16\t0.026845637583892617\t50.0\tvariant_caller\t-\t-\t-\tNo\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNC_000002.11\t29445216\t29445307\tATCAGGGCTTCCATGAGGAAATCCAGTTCGTCCTGTTCAGAGCACACTTCAGGCAGCGTCTGGGCAGAGAAGGGGAGGGTGGGGAGGAGGAG\tGCTGGTGCGGTCTCA\tALK:NM_004304:exon22:c.3451_3509TGAGACCGCACCAGC")
            #self.assertEqual(result[4].rstrip(),
            #    "sample1\tERBB2\tnonsynonymous_SNV\texon25\tp.G1015E\tc.G3044A\tNM_004448\t-\t2-other\tyes\tok\t1288\t3797\t613\t1971\t0.5190940215959968\t90.0\tvariant_caller\trs140272156\t0.004000000189989805\t0.0020000000949949026\tNo\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNC_000017.10\t37883141\t37883141\tG\tA\tERBB2:NM_001289937:exon25:c.G3044A:p.G1015E,ERBB2:NM_004448:exon25:c.G3044A:p.G1015E,ERBB2:NM_001005862:exon28:c.G2954A:p.G985E,ERBB2:NM_001289936:exon29:c.G2999A:p.G1000E")
            self.assertEqual(result[4].rstrip(),
                "sample1\tERBB2\tnonsynonymous_SNV\texon25\tp.G1015E\tc.G3044A\tNM_004448\t-\t4-other\tyes\tok\t1288\t3797\t1824\t1971\t0.5190940215959968\t90.0\tvariant_caller\trs140272156\t0.004000000189989805\t0.0020000000949949026\tNo\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNC_000017.10\t37883141\t37883141\tG\tA\tERBB2:NM_001289937:exon25:c.G3044A:p.G1015E,ERBB2:NM_004448:exon25:c.G3044A:p.G1015E,ERBB2:NM_001005862:exon28:c.G2954A:p.G985E,ERBB2:NM_001289936:exon29:c.G2999A:p.G1000E")

        report_2 = os.path.join(self.tempdir,"filtered.report_2")
        generate_filtered_mutations("sample1", "variant_caller", report_2, levels, self.hotspot_2, self.vcf_2 + ".gz" , self.pileup_2, self.reference, self.multibp, {"ERBB2": "NM_004448", "MET": "NM_001127500"})
        with open(report_2, 'r') as report_result:
            result = report_result.readlines()
            self.assertEqual(len(result),2)
            # Make sure that a snv isn't added to a indel hotspot
            self.assertEqual(result[1],
                "sample1\tKRAS\tstopgain\texon3\tp.G48X\tc.G142T\tNM_004985\t-\t4-other\tyes\tok\t4032\t4032\t4024\t8\t0.001984126984126984\t54.0\tvariant_caller\t-\t-\t-\tNo\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNC_000012.11\t25380316\t25380316\tC\tA\tKRAS:NM_004985:exon3:c.G142T:p.G48X,KRAS:NM_033360:exon3:c.G142T:p.G48X")

if __name__ == '__main__':
    import logging
    import sys
    logging.basicConfig(level=logging.CRITICAL, stream=sys.stdout, format='%(message)s')
    unittest.main()

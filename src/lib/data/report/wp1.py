from lib.data.files.hotspot import HotspotReader, ReportClass
from lib.data.files.models import ChrTranslater

from pysam import VariantFile, VCF

def generate_filtered_mutations(hotspot_file, vcf_file, pileup, chr_mapping):
    hotspot_reader = HotspotReader(hotspot_file)
    hotspots = [hotspot for hotspot in hotspot_reader]
    chr_translater = ChrTranslater(chr_mapping)

    variants = VCF()
    variants.connect(args.vcf)
    for variant in variants.

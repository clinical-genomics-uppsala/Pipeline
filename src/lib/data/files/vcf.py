# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from pysam import VariantFile

def add_AD_field_using_DP4(input_name, output_name):
    input_vcf = VariantFile(input_name,'r')
    input_vcf.header.info.add("AD", number=".", type='String', description="Allelic depths for the ref and alt alleles in the order listed")
    output_vcf = VariantFile(output_name,'w', header=input_vcf.header)
    #output_vcf.header.info.add("AD", number=".", type='Integer', description="Allelic depths for the ref and alt alleles in the order listed")
    for record in input_vcf:
        ref_pos, ref_neg, var_pos, var_neg = record.info['DP4']
        new_record = record.copy()
        new_record.info["AD"] = "{},{}".format(ref_pos + ref_neg, var_pos + var_neg)
        output_vcf.write(new_record)

def add_contigs_to_header(input_name, output_name, contig_file, assembly):
    from src.lib.data.files.reference import InfoImporter
    info = InfoImporter(contig_file)
    input_vcf = VariantFile(input_name,'r')
    for key in info:
            input_vcf.header.contigs.add(key, length=info[key]['length'], assembly=assembly)
    output_vcf = VariantFile(output_name,'w', header=input_vcf.header)
    #output_vcf.header.info.add("AD", number=".", type='Integer', description="Allelic depths for the ref and alt alleles in the order listed")
    for record in input_vcf:
        output_vcf.write(record)

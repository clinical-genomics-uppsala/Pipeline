import pysam
import src.lib.data.files.vcf as vcf

def get_annoation_data(variant,info_name):
    if isinstance(variant, pysam.VariantRecord):
        try:
            data = variant.info[info_name]
            if data is None:
                return "-"
            if isinstance(data,tuple):
                return ",".join(data)
            return data
        except KeyError:
            return "-"


def get_sample_value(variant, sample):
    if isinstance(variant, pysam.VariantRecord):
        return variant.samples['SAMPLE']


def get_vaf(variant,sample):
    if isinstance(variant, pysam.VariantRecord):
        #print(variant.formats['AD'])

        depth, ref, var = get_depth(variant, sample)
        return int(var)/(int(depth))

def get_depth(variant,sample):
    if isinstance(variant, pysam.VariantRecord):
        return (get_sample_value(variant,sample)['DP'],) + get_sample_value(variant,sample)['AD']

def get_report_type(variant):
    if isinstance(variant, pysam.VariantRecord):
        if vcf.is_indel(variant) or vcf.is_multibp_sub(variant):
            return "2-indel"
    return "1-hotspot"

def get_start_position(variant):
    if isinstance(variant, pysam.VariantRecord):
        if vcf.is_indel(variant):
            return variant.start + 2
        else:
            return variant.start + 1

def get_ref_sequence(variant):
    if isinstance(variant, pysam.VariantRecord):
        if len(variant.alleles) != 2:
            raise Exception("Unhandled exception")
        if vcf.is_indel(variant):
            if len(variant.alleles[0]) == 1:
                return "-"
            else:
                return variant.alleles[0][1:]
        else:
            return variant.alleles[0]

def get_var_sequence(variant):
    if isinstance(variant, pysam.VariantRecord):
        if len(variant.alleles) != 2:
            raise Exception("Unhandled exception")
        if vcf.is_indel(variant):
            if len(variant.alleles[1]) == 1:
                return "-"
            else:
                return variant.alleles[1][1:]
        else:
            return variant.alleles[1]

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from pysam import VariantFile

def check_overlapp(variant, chrom, position):
    return  variant.chrom == chrom and variant.start <= position <= variant.stop

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

def create_sample_format_from_info_lofreq(sample, input_name, output_name, skip_gt=False):
    input_vcf = VariantFile(input_name,'r')
    input_vcf.header.formats.add("AF", number=1, type='Float', description="Allele Frequency")
    input_vcf.header.formats.add("AD", number=".", type='String', description="Allelic sample depths for the ref and alt alleles in the order listed")
    input_vcf.header.formats.add("DP", number=1, type='Integer', description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)")
    input_vcf.header.formats.add("DP4", number=4, type='Integer', description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases")
    input_vcf.header.formats.add("GT", number=".", type="String", description="Genotype")

    input_vcf.header.add_sample(sample)
    output_vcf = VariantFile(output_name,'w', header=input_vcf.header)
    for record in input_vcf:
        ad = record.info["AD"]
        af = record.info["AF"]
        dp = 0
        for d in ad:
            dp = dp + int(d)
        fields = {"AF": af, "DP4": record.info["DP4"], "DP": dp, "AD": ad, "GT": (record.alleles[1],record.alleles[0])}
        new_record = output_vcf.new_record(record.chrom,record.start,record.stop,record.alleles, record.id, record.qual,record.filter,record.info,[fields]) #,
        output_vcf.write(new_record)


def process_annovar_refseq(data, gene, prefered_transcripts):
    genes = []
    transcripts = []
    exons_introns = []
    cds_changes = []
    aa_changes = []
    comm = "-"
    if isinstance(gene,list):
        if len(gene) != 1:
            raise Exception("Multiple genes: " + str(gene))
        else:
            gene = gene[0]
    else:
        if "," in gene:
            raise Exception("Multiple genes in variable: " + gene)
    if data and not data == ('-',):
        if gene in prefered_transcripts:
            found = False
            for d in data.split(","):
                if prefered_transcripts[gene] in d:
                    columns = d.split(":")
                    gene, transcript, exon_intron, cds_change, aa_change = "-", "-", "-", "-", "-"
                    if not d == ".":
                        if len(columns) == 5:
                            gene, transcript, exon_intron, cds_change, aa_change = columns
                        elif len(columns) == 4:
                            gene, transcript, exon_intron, cds_change = columns
                        else:
                            raise Exception("Unhandled number of data points when parsing refseq: {}".format(d))
                        found = True
                    genes.append(gene)
                    transcripts.append(transcript)
                    exons_introns.append(exon_intron)
                    cds_changes.append(cds_change)
                    aa_changes.append(aa_change)
            if not found:
                comm = "altTranscript"
        else:
            d = data.split(",")[0] #Why only pick the first?
            columns = d.split(":")
            gene, transcript, exon_intron, cds_change, aa_change = "-", "-", "-", "-", "-"
            if not d == ".":
                if len(columns) == 5:
                    gene, transcript, exon_intron, cds_change, aa_change = columns
                elif len(columns) == 4:
                    gene, transcript, exon_intron, cds_change = columns
                else:
                    raise Exception("Unhandled number of data points when parsing refseq: {}".format(d))
            genes.append(gene)
            transcripts.append(transcript)
            exons_introns.append(exon_intron)
            cds_changes.append(cds_change)
            aa_changes.append(aa_change)
    else:
        genes.append("-")
        transcripts.append("-")
        exons_introns.append("-")
        cds_changes.append("-")
        aa_changes.append("-")
    return (",".join(set(genes)), ",".join(set(transcripts)), ",".join(set(exons_introns)), ",".join(set(cds_changes)), ",".join(set(aa_changes)), comm)

def is_multibp_sub(variant):
    if len(variant.alleles) != 2:
        raise Exception("Unhandled case: " + str(variant.alleles))
    return len(variant.alleles[0]) > 1 and len(variant.alleles[1]) > 1


def is_indel(variant):
    if len(variant.alleles) != 2:
        raise Exception("Unhandled case: " + str(variant.alleles))
    return (len(variant.alleles[0]) == 1 and len(variant.alleles[1]) > 1) or (len(variant.alleles[1]) == 1 and len(variant.alleles[0]) > 1)

#def process_clinvar_data(data):
    #"""
    #CLINSIG=pathogenic,pathogenic|pathogenic|pathogenic|pathogenic|pathogenic|pathogenic|pathogenic
    #"""
    #clinvar = {}
    #for d in data.split("|"):
        #try:
            #key, value = d.split("=")
            #clinvar[key] = value
        #except ValueError:
            #continue
    #return clinvar

from src.lib.data.files.hotspot import Reader as HotspotReader
from src.lib.data.files.pileup import Reader as PileupReader
from src.lib.data.files.hotspot import ReportClass
from src.lib.data.files.models import ChrTranslater, Hotspot, MultiBpVariantData
import src.lib.data.files.utils as utils
import src.lib.data.files.vcf as vcf
from pysam import VariantFile, VCF, tabix_index

def _create_filtered_mutations_header():
    return "\t".join([
        "Sample",
        "Gene",
        "Variant_type",
        "Exon",
        "AA_change",
        "CDS_change",
        "Accession_number",
        "Comment",
        "Report",
        "Found",
        "Min_read_depth300",
        "Pileup_depth",
        "Total_read_depth",
        "Reference_read_depth",
        "Variant_read_depth",
        "Variant_allele_ratio",
        "dbSNP_id",
        "Ratio_1000G",
        "Ratio_ESP6500",
        "Clinically_flagged_dbSNP",
        "Cosmic",
        "ClinVar_CLNDB",
        "Clinval_CLINSIG",
        "Reference_plus_amplicons",
        "Reference_minus_amplicons",
        "Variant_plus_amplicons",
        "Variant_minus_amplicons",
        "Strands_A_F+F-S+S-",
        "Strands_G_F+F-S+S-",
        "Strands_C_F+F-S+S-",
        "Strands_T_F+F-S+S-",
        "Strands_Ins",
        "Strands_Del",
        "Ref_aligned_amplicons",
        "Var_aligned_amplicons",
        "Chr",
        "Start",
        "End",
        "Reference_base",
        "Variant_base",
        "All_transcripts_annotation"
    ])


def get_report_type(ref, var, insertion, deletion, mut_type):
    if ref == '-' or var == '-':
        return "2-indel"
    elif "+" in insertion or "+" in deletion:
        return "2-indel"
    return mut_type

def get_read_level(read_levels, rd):
    try:
        for (level, depth_status, analyzable) in read_levels:
            rd = int(rd)
            if rd >= int(level):
                return depth_status, analyzable
    except:
        pass
    return "-", "zero"

def _print_variant(writer, variant, sample, chrom, mutation_type, levels, pileup_depth, multibp, prefered_transcripts):
    depth_staus, found = get_read_level(levels, pileup_depth)
    variant_type = utils.get_annoation_data(variant, "Func.refGene")
    exonic_type = utils.get_annoation_data(variant, "ExonicFunc.refGene")
    genes = utils.get_annoation_data(variant, "Gene.refGene")
    if "splicing" in variant_type:
        exonic_type = variant_type
    _, transcripts, exons_introns, cds_changes, aa_changes, comm = vcf.process_annovar_refseq(utils.get_annoation_data(variant, "AAChange.refGene"), genes, prefered_transcripts)

    genes = utils.get_annoation_data(variant, "Gene.refGene")
    total_depth, ref_depth, var_depth = utils.get_depth(variant, sample)
    vaf = utils.get_vaf(variant, sample)
    #gene_details = utils.get_annoation_data(variant, "GeneDetail.refGene")
    dbSNP_id =  utils.get_annoation_data(variant, "snp138")
    nonflagged = utils.get_annoation_data(variant, "snp138NonFlagged")
    ratio_1000g = utils.get_annoation_data(variant, "1000g2015aug_eur")
    ratio_esp6500 = utils.get_annoation_data(variant, "esp6500siv2_ea")
    clinical_flagged = "Yes" if dbSNP_id.startswith("rs") and ("-" == nonflagged or nonflagged == ".") else "No"
    cosmic = utils.get_annoation_data(variant, "cosmic70")
    clinsig = utils.get_annoation_data(variant, "clinvar_20150629").replace("CLINSIG=","")
    clndb = utils.get_annoation_data(variant, "CLNDBN");
    #ToDo make sure vcf parse string correctly
    #ToDo, work around to handle incorrect pysam splitting  of info field
    occurrence = utils.get_annoation_data(variant, "OCCURENCE")
    if cosmic == ".":
        cosmic="-"
    else:
        cosmic = cosmic + ";" + "OCCURENCE=" + occurrence
    if clinsig == ".":
        clinsig = "-"
    #print(variant.info.items())
    if dbSNP_id == ".":
        dbSNP_id = "-"
    position = utils.get_start_position(variant)
    ref_seq = utils.get_ref_sequence(variant)
    var_seq = utils.get_var_sequence(variant)

    multibp_variant = multibp.get_data(chrom, position, variant.stop, ref_seq, var_seq)
    if multibp_variant is not None:
        aa_changes = multibp_variant.AA_CHANGE
        cds_changes = multibp_variant.CDS_CHANGE
        if not transcripts == multibp_variant.TRANSCRIPT:  # if there is a different accession number set comment to altTranscript
            transcripts = multibp_variant.TRANSCRIPT
            comm = "altTranscript"

    writer.write("\n" + (_construct_entry(sample=sample , gene=genes, variant_type=exonic_type, exon=exons_introns,
          aa_change=aa_changes, cds_change=cds_changes,
          accession_number = transcripts, comment=comm,
          report=mutation_type, found=found, min_read_depth300=depth_staus, pileup_depth=pileup_depth, total_read_depth=total_depth,
          reference_read_depth=ref_depth, variant_read_depth=var_depth, variant_allele_ratio=vaf,
          dbsnp_id=dbSNP_id, ratio_1000g=ratio_1000g, ratio_esp6500=ratio_esp6500, clinically_flagged_dbsnp=clinical_flagged,
          cosmic=cosmic, clinvar_clndb=clndb, clinval_clinsig=clinsig, reference_plus_amplicons="-",
          reference_minus_amplicons="-", variant_plus_amplicons="-" , variant_minus_amplicons="-",
          strands_a_ffss="-", strands_g_ffss="-", strands_c_ffss="-", strands_t_ffss="-", strands_ins="-",
          strands_del="-", ref_aligned_amplicons="-", var_aligned_amplicons="-",
          chr = chrom, start = position, end=variant.stop,
          reference_base=ref_seq, variant_base=var_seq, all_transcripts_annotation=utils.get_annoation_data(variant, "AAChange.refGene"))))

def _construct_entry(sample , gene, variant_type, exon, aa_change, cds_change, accession_number,
            comment, report, found, min_read_depth300, pileup_depth, total_read_depth, reference_read_depth,
            variant_read_depth, variant_allele_ratio,
            dbsnp_id, ratio_1000g, ratio_esp6500, clinically_flagged_dbsnp,
            cosmic, clinvar_clndb, clinval_clinsig, reference_plus_amplicons,
            reference_minus_amplicons, variant_plus_amplicons , variant_minus_amplicons,
            strands_a_ffss, strands_g_ffss, strands_c_ffss, strands_t_ffss, strands_ins, strands_del,
            ref_aligned_amplicons, var_aligned_amplicons,
            chr, start, end, reference_base, variant_base, all_transcripts_annotation):
  return ("{sample}\t{gene}\t{variant_type}\t{exon}\t{aa_change}\t{cds_change}" +
         "\t{accession_number}\t{comment}\t{report}\t{found}" +
         "\t{min_read_depth300}\t{pileup_depth}\t{total_read_depth}\t{reference_read_depth}" +
         "\t{variant_read_depth}\t{variant_allele_ratio}\t{dbsnp_id}" +
         "\t{ratio_1000g}\t{ratio_esp6500}\t{clinically_flagged_dbsnp}" +
         "\t{cosmic}\t{clinvar_clndb}\t{clinval_clinsig}" +
         "\t{reference_plus_amplicons}\t{reference_minus_amplicons}" +
         "\t{variant_plus_amplicons}\t{variant_minus_amplicons}" +
         "\t{strands_a_ffss}\t{strands_g_ffss}\t{strands_c_ffss}" +
         "\t{strands_t_ffss}\t{strands_ins}\t{strands_del}" +
         "\t{ref_aligned_amplicons}\t{var_aligned_amplicons}" +
         "\t{chr}\t{start}\t{end}\t{reference_base}\t{variant_base}\t{all_transcripts_annotation}").format(
          sample=sample,
          gene = gene,
          variant_type = variant_type,
          exon = exon,
          aa_change = aa_change,
          cds_change = cds_change,
          accession_number = accession_number,
          comment = comment,
          report = report,
          found = found,
          min_read_depth300 = min_read_depth300,
          pileup_depth = pileup_depth,
          total_read_depth = total_read_depth,
          reference_read_depth = reference_read_depth,
          variant_read_depth = variant_read_depth,
          variant_allele_ratio = variant_allele_ratio,
          dbsnp_id = dbsnp_id,
          ratio_1000g = ratio_1000g,
          ratio_esp6500 = ratio_esp6500,
          clinically_flagged_dbsnp = clinically_flagged_dbsnp,
          cosmic = cosmic,
          clinvar_clndb = clinvar_clndb,
          clinval_clinsig = clinval_clinsig,
          reference_plus_amplicons = reference_plus_amplicons,
          reference_minus_amplicons = reference_minus_amplicons,
          variant_plus_amplicons = variant_plus_amplicons,
          variant_minus_amplicons = variant_minus_amplicons,
          strands_a_ffss = strands_a_ffss,
          strands_g_ffss = strands_g_ffss,
          strands_c_ffss = strands_c_ffss,
          strands_t_ffss = strands_t_ffss,
          strands_ins = strands_ins,
          strands_del = strands_del,
          ref_aligned_amplicons = ref_aligned_amplicons,
          var_aligned_amplicons = var_aligned_amplicons,
          chr = chr,
          start = start,
          end = end,
          reference_base = reference_base,
          variant_base = variant_base,
          all_transcripts_annotation = all_transcripts_annotation)

def _print_report(writer, sample, hotspot, other, levels, chr_translater, multibp, prefered_transcripts):
    for index, depth_variant in enumerate(hotspot.DEPTH_VARIANTS):
        #print(str(index) + " " + str(depth_variant))
        if not depth_variant['variants'] and not depth_variant['extended']:
            depth_staus, found = get_read_level(levels, depth_variant['depth'])
            if found == "yes":
                found = "no"
            writer.write("\n" + (_construct_entry(sample=sample , gene=hotspot.GENE, variant_type="-", exon=hotspot.EXON,
                      aa_change=hotspot.AA_MUTATION_SYNTAX, cds_change=hotspot.CDS_MUTATION_SYNTAX,
                      accession_number = hotspot.ACCESSION_NUMBER, comment=hotspot.COMMENT,
                      report=hotspot.REPORT, found=found, min_read_depth300=depth_staus, pileup_depth=depth_variant['depth'], total_read_depth="-",
                      reference_read_depth="-", variant_read_depth="-", variant_allele_ratio="-",
                      dbsnp_id="-", ratio_1000g="-", ratio_esp6500="-", clinically_flagged_dbsnp="-",
                      cosmic="-", clinvar_clndb="-", clinval_clinsig="-", reference_plus_amplicons="-",
                      reference_minus_amplicons="-", variant_plus_amplicons="-" , variant_minus_amplicons="-",
                      strands_a_ffss="-", strands_g_ffss="-", strands_c_ffss="-", strands_t_ffss="-", strands_ins="-",
                      strands_del="-", ref_aligned_amplicons="-", var_aligned_amplicons="-",
                      chr = hotspot.CHROMOSOME, start=(hotspot.START+index), end=(hotspot.START+index),
                      reference_base="-", variant_base="-", all_transcripts_annotation="-")))
        else:
            for var in depth_variant['variants']:
                #print("Depth: " + str(depth_variant['depth']) + "\t" + str(var))
                _print_variant(writer, var, sample, chr_translater.get_nc_value(var.chrom), utils.get_report_type(var), levels, depth_variant['depth'], multibp, prefered_transcripts)


def generate_filtered_mutations(sample, report, levels, hotspot_file, vcf_file, pileup_file, chr_mapping, multibp_file, prefered_transcripts):
    hotspot_reader = HotspotReader(hotspot_file)
    hotspots = [hotspot for hotspot in iter(hotspot_reader)]
    multibp = MultiBpVariantData(multibp_file)
    chr_translater = ChrTranslater(chr_mapping)
    pileup_reader = PileupReader(pileup_file)
    variants = VariantFile(vcf_file)
    other = []
    for variant in variants:
        if not len(variant.alts) == 1:
            raise Exception("Multiple allele found: " + str(variant.alts))
        added = False
        for hotspot in hotspots:
            if hotspot.add_variant(variant, chr_translater):
                added = True
        if not added:
            other.append((variant,"-"))
    for data in pileup_reader:
        added = False
        for hotspot in hotspots:
            if hotspot.add_depth(data, chr_translater):
                added = True
        for i in range(0,len(other)):
            #print(data.CHROMOSOME + "\t" + data.POSITION + "\n")
            if vcf.check_overlapp(other[i][0], data.CHROMOSOME, data.POSITION):
                other[i] = (other[i][0], data.DEPTH)

    with open(report,"w") as filtered_mutations:
        filtered_mutations.write(_create_filtered_mutations_header())
        #print("Hotspot")
        for hotspot in hotspots:
            _print_report(filtered_mutations, sample, hotspot, other, levels, chr_translater, multibp, prefered_transcripts)
        #print("Other")
        for var in other:
            _print_variant(filtered_mutations, var[0], sample, chr_translater.get_nc_value(var[0].chrom), "4-other", levels, var[1], multibp, prefered_transcripts)

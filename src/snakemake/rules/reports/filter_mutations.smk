# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

_generate_filtered_mutations_wp1_input_vcf = "variants/{sample}.{part}.{caller}.annovar.vcf"
try:
    _generate_filtered_mutations_wp1_input_vcf = generate_filtered_mutations_wp1_input_vcf
except:
    pass

_generate_filtered_mutations_wp1_input_pileup = "variants/{sample}.{part}.{caller}.annovar.pileup"
try:
    _generate_filtered_mutations_wp1_input_pileup = generate_filtered_mutations_wp1_input_pileup
except:
    pass

_generate_filtered_mutations_wp1_output = "reports/{sample}.{part}.{caller}.filteredMutation.tsv"
try:
    _generate_filtered_mutations_wp1_output = generate_filtered_mutations_wp1_output
except:
    pass

from src.lib.data.report.wp1 import generate_filtered_mutations

rule generate_filtered_mutations_wp1:
  input:
      vcf=_generate_filtered_mutations_wp1_input_vcf,
      pileup=_generate_filtered_mutations_wp1_input_pileup
  output:
      _generate_filtered_mutations_wp1_output
  log:
      "logs/filter_mutations/{sample}.{part}.{caller}.log"
  params:
      levels = config['depth_levels'],
      hotspot = lambda wildcards:  samples['hotspot'][wildcards.sample],
      chr_mapping = config['chr_mapping_file'],
      prefered_transcripts = config['chr_mapping_file'],
      multibp_file = config['multibp_file']
  run:
    generate_filtered_mutations(
        wildcards.sample,
        wildcards.caller,
        output[0],
        params.levels,
        params.hotspot,
        input.vcf,
        input.pileup,
        params.chr_mapping,
        params.multibp_file,
        params.prefered_transcripts)

import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

def generate_file_output_bam():
    return [os.path.join("alignment", str(row.Index)+ ending) for row in samples.itertuples() for ending in [".0000-consensus3.primerclip.bam"]]


def generate_file_output_fastq():
    return [os.path.join("fastq", str(row.Index)+ ending) for row in samples.itertuples() for ending in [".0000.consensus.1.R1.fq", ".0000.consensus.1.R2.fq"]]

def generate_file_output_vcf():
    return [os.path.join("variants", str(row.Index)+ ending) for row in samples.itertuples() for ending in [".0000-consensus3.lofreq.vcf"]]

rule all:
    input:
        generate_file_output_vcf() #generate_file_output_fastq() + generate_file_output_bam()
        #[*generate_file_output_fastq(), *generate_file_output_bam()]
        #["alignment/JI-2.L001-0003.bam", "alignment/JI-2.L001-0004.bam", "alignment/JI-2.L001-0000.bam", "alignment/JI-2.L001-0001.bam", "alignment/JI-2.L001-0002.bam"] #generate_file_output_bam()

include: "src/snakemake/workflow/swift.hs.snakemake"

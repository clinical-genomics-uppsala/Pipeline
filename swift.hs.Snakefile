import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

qc_files = [("fastqc", "_fastqc.html"),("fastqc" ,"_fastqc.zip")]
samstats_files = [("samstats", "0000.bam.txt"), ("samstats", "0000-consensus3.primerclip.rg.bam.txt")]

def generate_file_output_bam():
    return [os.path.join("alignment", str(row.Index)+ ending) for row in samples.itertuples() for ending in [".0000-consensus3.primerclip.bam"]]


def generate_file_output_fastq():
    return [os.path.join("fastq", str(row.Index)+ ending) for row in samples.itertuples() for ending in [".0000.consensus.1.R1.fq", ".0000.consensus.1.R2.fq"]]

def generate_file_output_vcf():
    return [os.path.join("variants", str(row.Index)+ ending) for row in samples.itertuples() for ending in [".0000-consensus3.lofreq.annovar.vcf"]]

def generate_filtered_mutations():
    return [os.path.join("reports", str(row.Index)+ ".0000-consensus3.lofreq.filteredMutation.tsv") for row in samples.itertuples()]

def generate_qc_reports():
    return [os.path.join("qc", f[0],  str(row.Index[0]) + "." + str(row.Index[1]) + "." + f[1])
        for row in units.itertuples()
            for f in qc_files]

def generate_samstats_reports():
    return [os.path.join("qc", f[0],  str(row.Index[0]) + "." + f[1])
        for row in units.itertuples()
            for f in samstats_files]

rule all:
    input:
        generate_filtered_mutations() + ["reports/multiqc.html"] #generate_file_output_fastq() + generate_file_output_bam()
        #[*generate_file_output_fastq(), *generate_file_output_bam()]
        #["alignment/JI-2.L001-0003.bam", "alignment/JI-2.L001-0004.bam", "alignment/JI-2.L001-0000.bam", "alignment/JI-2.L001-0001.bam", "alignment/JI-2.L001-0002.bam"] #generate_file_output_bam()

include: "src/snakemake/workflow/swift.hs.snakemake"

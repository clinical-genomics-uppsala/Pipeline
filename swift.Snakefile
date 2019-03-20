import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

def generate_filtered_mutations():
    return [os.path.join("reports", str(row.Index)+ ".0000.lofreq.filteredMutation.tsv") for row in samples.itertuples()]

rule all:
    input:
        generate_filtered_mutations() + ["reports/multiqc.html"] #generate_file_output_fastq() + generate_file_output_bam()
        #[*generate_file_output_fastq(), *generate_file_output_bam()]
        #["alignment/JI-2.L001-0003.bam", "alignment/JI-2.L001-0004.bam", "alignment/JI-2.L001-0000.bam", "alignment/JI-2.L001-0001.bam", "alignment/JI-2.L001-0002.bam"] #generate_file_output_bam()

include: "src/snakemake/workflow/swift.snakemake"

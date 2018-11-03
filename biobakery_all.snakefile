# set options with "--config "

import os

configfile: "config.yaml"

input_folder = config["input_folder"]
output_folder = config["output_folder"]


SAMPLES, = glob_wildcards(os.path.join(input_folder, "{samples}_R1_001.fastq.gz"))

include: "workflows/kneaddata.snakefile"
include: "workflows/metaphlan2.snakefile"
include: "workflows/humann2.snakefile"


rule all:
    input:
        os.path.join(output_folder, "report.html")


rule report:
    input:
        kneaddata = os.path.join(output_folder, "kneaddata/kneaddata_report.html"),
        metaphlan2 = os.path.join(output_folder, "metaphlan2/metaphlan2_report.html"),
        humann2 = os.path.join(output_folder, "humann2/humann2_report.html")
    output:
        os.path.join(output_folder, "report.html")
    run:
        from snakemake.utils import report
        report("""
        Pipeline works!!!
        """, output[0])

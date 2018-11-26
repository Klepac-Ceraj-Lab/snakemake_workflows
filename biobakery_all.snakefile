# set options with "--config "

import os

configfile: "config.yaml"

input_folder = config["input_folder"]
output_folder = config["output_folder"]

(samples, lanes) = glob_wildcards(os.path.join(input_folder, "{sample}_{lane,L\d+}_R1_001.fastq.gz"))
samples = list(set(samples))
samples.sort()
lanes = list(set(lanes))
lanes.sort()

include: "workflows/kneaddata.snakefile"
include: "workflows/metaphlan2.snakefile"
include: "workflows/humann2.snakefile"

kneadfolder = os.path.join(output_folder, "kneaddata")
metaphlanfolder = os.path.join(output_folder, "metaphlan2")
humannfolder = os.path.join(input_folder, "humann2")

rule all:
    input:
        os.path.join(output_folder, "report.html")


rule report:
    input:
        kneaddata = os.path.join(kneadfolder, "kneaddata_report.html"),
        metaphlan2 = os.path.join(metaphlanfolder, "metaphlan2_report.html"),
        humann2 = os.path.join(humannfolder, "humann2_report.html")
    output:
        os.path.join(output_folder, "report.html")
    run:
        from snakemake.utils import report
        report("""
        Pipeline works!!!
        """, output[0])

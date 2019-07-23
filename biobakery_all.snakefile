# set options with "--config "

import os, glob

configfile: "config.yaml"

input_folder = config["input_folder"]
output_folder = config["output_folder"]
log_folder = os.path.join(output_folder, "logs")


if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

if not os.path.isdir(log_folder):
    os.mkdir(log_folder)

(samples, lanes) = glob_wildcards(os.path.join(input_folder, "{sample}_{lane,L\d+}_R1_001.fastq.gz"))
samples = list(set(samples))
samples.sort()
lanes = list(set(lanes))
lanes.sort()

print(samples)

kneadfolder = os.path.join(output_folder, "kneaddata")
metaphlanfolder = os.path.join(output_folder, "metaphlan2")
humannfolder = os.path.join(output_folder, "humann2")

include: "workflows/kneaddata.snakefile"
include: "workflows/metaphlan2.snakefile"
include: "workflows/humann2.snakefile"

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

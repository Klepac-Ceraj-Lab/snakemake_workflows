import os, glob

input_folder = config["input_folder"]
output_folder = config["output_folder"]
log_folder = os.path.join(output_folder, "logs")

if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

if not os.path.isdir(log_folder):
    os.mkdir(log_folder)

kneadfolder = os.path.join(output_folder, "kneaddata")
metaphlanfolder = os.path.join(output_folder, "metaphlan2")
humannfolder = os.path.join(output_folder, "humann2")

samples = glob_wildcards(os.path.join(metaphlanfolder, "main", "{sample}_profile.tsv"))
samples = list(set(samples))
samples.sort()

include: "workflows/humann2.snakefile"

rule all:
    input:
        os.path.join(humannfolder, "humann2_report.html")

# set options with "--config "

import os, glob

configfile: "config.yaml"

include: "setup/directories.snakefile"
include: "setup/echo_samples.snakefile"


kneadfolder = os.path.join(output_folder, "kneaddata")
include: "workflows/kneaddata.snakefile"

metaphlanfolder = os.path.join(output_folder, "metaphlan")
include: "workflows/metaphlan.snakefile"

rule all:
    input:
        os.path.join(output_folder, "report.html")


rule report:
    input:
        kneaddata = os.path.join(kneadfolder, "kneaddata_report.html"),
        metaphlan = os.path.join(metaphlanfolder, "metaphlan_report.html"),
    output:
        os.path.join(output_folder, "report.html")
    run:
        from snakemake.utils import report
        report("""
        Pipeline works!!!
        """, output[0])

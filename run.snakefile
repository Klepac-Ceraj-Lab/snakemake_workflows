# set options with "--config "
import os, glob

configfile: "config.yaml"

include: "setup/directories.snakefile"
include: config["samples"]

if config["run"] == "all":
    include: "run_all.snakefile"
elif config["run"] == "kneaddata":
    include: "run_kneaddata.snakefile"
elif config["run"] == "metaphlan":
    include: "run_metaphlan.snakefile"
else:
    raise ValueError("Unknown workflow: {}".format(config["run"]))


rule all:
    input:
        os.path.join(output_folder, "report.html")


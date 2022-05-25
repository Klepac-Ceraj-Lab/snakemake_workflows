# set options with "--config "
import os, glob

configfile: "config.yaml"

include: "setup/directories.smk"
include: config["samples"]

if config["run"] == "all":
    include: "run_all.smk"
elif config["run"] == "kneaddata":
    include: "run_kneaddata.smk"
elif config["run"] == "metaphlan":
    include: "run_metaphlan.smk"
else:
    raise ValueError("Unknown workflow: {}".format(config["run"]))


rule all:
    input:
        os.path.join(output_folder, "report.html")

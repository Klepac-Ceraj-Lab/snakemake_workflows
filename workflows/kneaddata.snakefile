kneadfolder = os.path.join(output_folder, "kneaddata")

if config["reads"] == "paired":
    include: "kneaddata_paired.snakefile"
elif config["reads"] == "single":
    include: "kneaddata_single.snakefile"


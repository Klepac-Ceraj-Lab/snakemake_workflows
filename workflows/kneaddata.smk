kneadfolder = os.path.join(output_folder, "kneaddata")

if config["reads"] == "paired":
    include: "kneaddata_paired.smk"
elif config["reads"] == "single":
    include: "kneaddata_single.smk"


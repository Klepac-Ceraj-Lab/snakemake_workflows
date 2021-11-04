rule kneaddata_cat_single:
    input: os.path.join(input_folder, "{sample}.fastq.gz")
    output: os.path.join(input_folder, "{sample}.fastq.gz")
    run:
        shell("cat {input} > {output}")

rule kneaddata:
    input:
        seq = os.path.join(input_folder, "{sample}.fastq.gz"),
        db = config["databases"]["human_sequences"]
    output:
        seq = temp(os.path.join(kneadfolder, "{sample}_kneaddata.fastq")),
        log = os.path.join(kneadfolder, "{sample}_kneaddata.log")
    run:
        shell("kneaddata --input {{input.seq}} --reference-db {{input.db}} --output {} --output-prefix {{wildcards.sample}}_kneaddata".format(kneadfolder))

rule kneaddata_counts:
    input: expand(os.path.join(kneadfolder, "{sample}_kneaddata.log"), sample = samples)
    output: os.path.join(kneadfolder, "kneaddata_read_counts.txt")
    shell:
        "kneaddata_read_count_table --input {} --output {{output}}".format(kneadfolder)

rule kneaddata_report:
    input:
        counts = os.path.join(kneadfolder, "kneaddata_read_counts.txt"),
        fwd = expand(os.path.join(kneadfolder, "{sample}_kneaddata.fastq"), sample = samples),
    output:
        os.path.join(kneadfolder, "kneaddata_report.html")
    run:
        from snakemake.utils import report
        report("""
        KneadData works!!!
        """, output[0])

kneadfolder = os.path.join(output_folder, "kneaddata/kneaddata_output/")
samples, = glob_wildcards(os.path.join(input_folder, "{samples}_R1_001.fastq.gz"))

rule kneadfiles:
    input:
        fwd = expand(os.path.join(kneadfolder, "{sample}_R1_001_kneaddata_paired_1.fastq"), sample = samples),
        rev = expand(os.path.join(kneadfolder, "{sample}_R1_001_kneaddata_paired_2.fastq"), sample = samples)


rule kneaddata_fwd:
    input:
        fwd = os.path.join(input_folder, "{samples}_R1_001.fastq.gz"),
        rev = os.path.join(input_folder, "{samples}_R2_001.fastq.gz"),
        db = config["databases"]["human_sequences"]
    output:
        fwd = os.path.join(kneadfolder, "{samples}_R1_001_kneaddata_paired_1.fastq"),
        rev = os.path.join(kneadfolder, "{samples}_R1_001_kneaddata_paired_2.fastq")
    run:
        shell("kneaddata --input {{input.fwd}} --input {{input.rev}} --reference-db {{input.db}} --output {}".format(kneadfolder))


rule kneaddata_counts:
    input:
        fwd = expand(os.path.join(kneadfolder, "{sample}_R1_001_kneaddata_paired_1.fastq"), sample = samples),
        rev = expand(os.path.join(kneadfolder, "{sample}_R1_001_kneaddata_paired_2.fastq"), sample = samples)
    output:
        os.path.join(output_folder, "kneaddata/kneaddata_read_counts.txt")
    shell:
        "kneaddata_read_count_table --input {} --output {{output}}".format(kneadfolder)


rule kneaddata_report:
    input:
        os.path.join(output_folder, "kneaddata/kneaddata_read_counts.txt")
    output:
        os.path.join(output_folder, "kneaddata/kneaddata_report.html")
    run:
        from snakemake.utils import report
        report("""
        KneadData works!!!
        """, output[0])

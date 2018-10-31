rule kneaddata_filter:
    input:
        fwd = expand(os.path.join(input_folder, "{samples}_R1_001.fastq.gz"), samples = SAMPLES),
        rev = expand(os.path.join(input_folder, "{samples}_R2_001.fastq.gz"), samples = SAMPLES),
        db = config["databases"]["human_sequences"]
    output:
        folder = os.path.join(output_folder, "kneaddata/kneaddata_output/"),
        fwd = expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_1.fastq"), samples = SAMPLES),
        rev = expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_2.fastq"), samples = SAMPLES)
    run:
        for f,r in zip(input.fwd,input.rev):
            shell("kneaddata --input {f} --input {r} --reference-db {input.db} --output {output.folder}")


rule kneaddata_counts:
    input:
        os.path.join(output_folder, "kneaddata/kneaddata_output/")
    output:
        os.path.join(output_folder, "kneaddata/kneaddata_read_counts.txt")
    shell:
        "kneaddata_read_count_table --input {input} --output {output}"


rule kneaddata_report:
    input:
        filter = os.path.join(output_folder, "kneaddata/kneaddata_output/"),
        counts = os.path.join(output_folder, "kneaddata/kneaddata_read_counts.txt")
    output:
        os.path.join(output_folder, "kneaddata/kneaddata_report.html")
    run:
        from snakemake.utils import report
        report("""
        KneadData works!!!
        """, output[0])
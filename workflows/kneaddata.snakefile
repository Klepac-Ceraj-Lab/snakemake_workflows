rule kneaddata_cat1:
    input: expand(os.path.join(input_folder, "{{sample}}_{lane}_R1_001.fastq.gz"), lane=lanes)
    output: os.path.join(input_folder, "{sample}_1.fastq.gz")
    run:
        shell("cat {input} > {output}")

rule kneaddata_cat2:
    input: expand(os.path.join(input_folder, "{{sample}}_{lane}_R2_001.fastq.gz"), lane=lanes)
    output: os.path.join(input_folder, "{sample}_2.fastq.gz")
    run:
        shell("cat {input} > {output}")

rule kneaddata:
    input:
        fwd = os.path.join(input_folder, "{sample}_1.fastq.gz"),
        rev = os.path.join(input_folder, "{sample}_2.fastq.gz"),
        db = config["databases"]["human_sequences"]
    output:
        fwd = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"),
        rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq")
    run:
        shell("kneaddata --input {{input.fwd}} --input {{input.rev}} --reference-db {{input.db}} --output {} --output-prefix {{wildcards.sample}}_kneaddata".format(kneadfolder))

rule kneaddata_cat_result:
    input:
        fwd = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"),
        rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq")
    output: os.path.join(kneadfolder, "{sample}.fastq")
    run:
        shell("cat {input} > {output}")

# rule kneaddata_gzip:
#     input:
#         fwd = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"),
#         rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq")
#     output:
#         fwd = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq.gz"),
#         rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq.gz")
#     run:
#          for f in glob.glob(os.path.join(kneadfolder, "{wildcards.sample}*")):
#              shell("gzip -v {}".format(f))
#
rule kneaddata_counts:
    input:
        fwd = expand(os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"), sample=samples),
        rev = expand(os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq"), sample=samples)
    output:
        os.path.join(kneadfolder, "kneaddata_read_counts.txt")
    shell:
        "kneaddata_read_count_table --input {} --output {{output}}".format(kneadfolder)

rule kneaddata_report:
    input:
        counts = os.path.join(kneadfolder, "kneaddata_read_counts.txt"),
        merged_fastqs = expand(os.path.join(kneadfolder, "{sample}.fastq"), sample=samples)
    output:
        os.path.join(kneadfolder, "kneaddata_report.html")
    run:
        from snakemake.utils import report
        report("""
        KneadData works!!!
        """, output[0])

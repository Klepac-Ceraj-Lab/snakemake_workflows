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
        fwd = temp(dynamic(expand(os.path.join(kneadfolder, "{sample}_{{filter_type}}_1.fastq"), sample = samples))),
        rev = temp(dynamic(expand(os.path.join(kneadfolder, "{sample}_{{filter_type}}_2.fastq"), sample = samples)))
    run:
        shell("kneaddata --input {{input.fwd}} --input {{input.rev}} --reference-db {{input.db}} --output {} --output-prefix {{wildcards.sample}}_kneaddata".format(kneadfolder))

rule kneaddata_gzip:
    input:
        fwd = dynamic(os.path.join(kneadfolder, "{sample}_kneaddata{{filter_type}}.fastq")),
        rev = dynamic(os.path.join(kneadfolder, "{sample}_kneaddata{{filter_type}}.fastq"))
    output:
        fwd = dynamic(os.path.join(kneadfolder, "{sample}_kneaddata{{filter_type}}.fastq.gz")),
        rev = dynamic(os.path.join(kneadfolder, "{sample}_kneaddata{{filter_type}}.fastq.gz"))
    run:
         shell("gzip {input.fwd}")
         shell("gzip {input.rev}")

rule kneaddata_counts:
    input:
        fwd = dynamic(expand(os.path.join(kneadfolder, "{sample}_{{filter_type}}_1.fastq"), sample = samples)),
        rev = dynamic(expand(os.path.join(kneadfolder, "{sample}_{{filter_type}}_2.fastq"), sample = samples))
    output:
        os.path.join(kneadfolder, "kneaddata_read_counts.txt")
    shell:
        "kneaddata_read_count_table --input {} --output {{output}}".format(kneadfolder)

rule kneaddata_gzip_foward:
    input: dynamic(expand(os.path.join(kneadfolder, "{sample}_{{filter_type}}_1.fastq"), sample = samples))
    output: dynamic(expand(os.path.join(kneadfolder, "{sample}_{{filter_type}}_1.fastq.gz"), sample = samples))
    run:
        shell("gzip {input}")

rule kneaddata_gzip_reverse:
    input: dynamic(expand(os.path.join(kneadfolder, "{sample}_{{filter_type}}_2.fastq"), sample = samples))
    output: dynamic(expand(os.path.join(kneadfolder, "{sample}_{{filter_type}}_2.fastq.gz"), sample = samples))
    run:
        shell("gzip {input}")

rule kneaddata_report:
    input:
        fwd: dynamic(expand(os.path.join(kneadfolder, "{sample}_{{filter_type}}_1.fastq.gz"), sample = samples))
        rev: dynamic(expand(os.path.join(kneadfolder, "{sample}_{{filter_type}}_2.fastq.gz"), sample = samples))
    output:
        os.path.join(kneadfolder, "kneaddata_report.html")
    run:
        from snakemake.utils import report
        report("""
        KneadData works!!!
        """, output[0])

rule kneaddata_cat_pair1:
    input: expand(os.path.join(input_folder, "{{sample}}_{lane}_R1_001.fastq.gz"), lane=lanes)
    output: os.path.join(input_folder, "{sample}_1.fastq.gz")
    run:
        shell("cat {input} > {output}")
rule kneaddata_cat_pair2:
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
        rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq"),
        log = os.path.join(kneadfolder, "{sample}_kneaddata.log")
    run:
        shell("kneaddata --input {{input.fwd}} --input {{input.rev}} --reference-db {{input.db}} --output {} --output-prefix {{wildcards.sample}}_kneaddata --trimmomatic /home/vklepacc/miniconda3/share/trimmomatic-0.39-2/".format(kneadfolder))

rule metaphlan_cat:
    input:
        fwd = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"),
        rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq"),
    output: temp(os.path.join(kneadfolder, "{sample}_kneaddata.fastq"))
    run:
        shell("cat {input} > {output}")


rule kneaddata_counts:
    input: expand(os.path.join(kneadfolder, "{sample}_kneaddata.log"), sample = samples)
    output: os.path.join(kneadfolder, "kneaddata_read_counts.txt")
    shell:
        "kneaddata_read_count_table --input {} --output {{output}}".format(kneadfolder)

rule kneaddata_report:
    input:
        counts = os.path.join(kneadfolder, "kneaddata_read_counts.txt"),
        fwd = expand(os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"), sample = samples),
        rev = expand(os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq"), sample = samples)
    output:
        os.path.join(kneadfolder, "kneaddata_report.html")
    run:
        from snakemake.utils import report
        report("""
        KneadData works!!!
        """, output[0])

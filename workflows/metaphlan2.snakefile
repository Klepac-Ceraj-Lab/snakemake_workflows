rule metaphlan2_cat:
    input:
        fwd = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"),
        rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq"),
    output: temp(os.path.join(kneadfolder, "{sample}.fastq"))
    run: shell("cat {input} > {output}")

rule metaphlan2_reads:
    input: os.path.join(kneadfolder, "{sample}.fastq")
    output:
        profile = os.path.join(output_folder, "metaphlan2/main/{sample}_profile.txt"),
        bowtie = os.path.join(output_folder, "metaphlan2/main/{sample}_bowtie2.txt"),
        sam = os.path.join(output_folder, "metaphlan2/main/{sample}.sam.bz2")
    run:
        shell("metaphlan2.py {input} {output.profile} --bowtie2out {output.bowtie} --samout {output.sam} --input_type fastq --nproc 8") # TODO: get nproc from settings


rule metaphlan2_merge:
    input:
        expand(os.path.join(output_folder, "metaphlan2/main/{sample}_profile.txt"), sample = SAMPLES)
    output:
        os.path.join(output_folder, "metaphlan2/merged/merged_abundance_table.txt")
    shell:
        "merge_metaphlan_tables.py {input} > {output}"


rule metaphlan2_report:
    input:
        os.path.join(output_folder, "metaphlan2/merged/merged_abundance_table.txt")
    output:
        os.path.join(output_folder, "metaphlan2/metaphlan2_report.html")
    run:
        from snakemake.utils import report
        report("""
        MetaPhlAn2 works!!!
        """, output[0])

rule metaphlan2_cat:
    input:
        fwd = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"),
        rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq"),
    output: temp(os.path.join(kneadfolder, "{sample}_merged.fastq"))
    run:
        shell("cat {input} > {output}")

rule metaphlan2:
    input: os.path.join(kneadfolder, "{sample}_merged.fastq")
    output:
        profile = os.path.join(metaphlanfolder, "main", "{sample}_profile.tsv"),
        bowtie = os.path.join(metaphlanfolder, "main", "{sample}_bowtie2.tsv"),
        sam = temp(os.path.join(metaphlanfolder, "main", "{sample}.sam"))
    run:
        shell("metaphlan2.py {input} {output.profile} --bowtie2out {output.bowtie} --samout {output.sam} --input_type fastq --nproc 8") # TODO: get nproc from settings


rule metaphlan2_merge:
    input:
        expand(os.path.join(metaphlanfolder, "main", "{sample}_profile.tsv"), sample = samples)
    output:
        os.path.join(metaphlanfolder, "merged", "merged_abundance_table.tsv")
    run:
        shell("merge_metaphlan_tables.py {input} > {output}")

rule metaphlan2_bz2:
    input: os.path.join(metaphlanfolder, "main", "{sample}.sam")
    output: os.path.join(metaphlanfolder, "main", "{sample}.sam.bz2")
    run:
        shell("bzip2 {input}")

rule metaphlan2_report:
    input:
        abundance_table = os.path.join(metaphlanfolder, "merged", "merged_abundance_table.tsv"),
        sams = expand(os.path.join(metaphlanfolder, "main", "{sample}.sam.bz2"), sample = samples)
    output:
        os.path.join(metaphlanfolder, "metaphlan2_report.html")
    run:
        from snakemake.utils import report
        report("""
        MetaPhlAn2 works!!!
        """, output[0])

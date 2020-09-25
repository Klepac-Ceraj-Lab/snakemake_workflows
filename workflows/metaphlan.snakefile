metaphlanfolder = os.path.join(output_folder, "metaphlan")

rule metaphlan:
    input: os.path.join(kneadfolder, "{sample}_kneaddata.fastq")
    output:
        profile = os.path.join(metaphlanfolder, "{sample}_profile.tsv"),
        bowtie = os.path.join(metaphlanfolder, "{sample}_bowtie2.tsv"),
        sam = temp(os.path.join(metaphlanfolder, "{sample}.sam"))
    run:
        shell("metaphlan {input} {output.profile} --bowtie2out {output.bowtie} --samout {output.sam} --input_type fastq --nproc 8 --bowtie2db /nobackup1/vklepacc/databases/chocophlan_markers") # TODO: get nproc from settings

# rule metaphlan_bz2:
#     input: os.path.join(metaphlanfolder, "{sample}.sam")
#     output: os.path.join(metaphlanfolder, "{sample}.sam.bz2")
#     run:
#         shell("bzip2 {input}")

rule metaphlan_report:
    input:
        profiles = expand(os.path.join(metaphlanfolder, "{sample}_profile.tsv"), sample = samples),
        sams = expand(os.path.join(metaphlanfolder, "{sample}.sam"), sample = samples)
    output:
        os.path.join(metaphlanfolder, "metaphlan_report.html")
    run:
        from snakemake.utils import report
        report("""
        metaphlan works!!!
        """, output[0])

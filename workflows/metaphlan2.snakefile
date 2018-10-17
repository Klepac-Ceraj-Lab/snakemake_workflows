rule metaphlan2_reads:
    input:
        fwd = expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_1.fastq"), samples = SAMPLES),
        rev = expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_2.fastq"), samples = SAMPLES)
    output:
        profiles = expand(os.path.join(output_folder, "metaphlan2/main/{samples}_profile.txt"), samples = SAMPLES),
        bowties = expand(os.path.join(output_folder, "metaphlan2/main/{samples}_bowtie2.txt"), samples = SAMPLES),
        sams = expand(os.path.join(output_folder, "metaphlan2/main/{samples}.sam.bz2"), samples = SAMPLES)
    run:
        for f,r,p,b,s in zip(input.fwd,input.rev,output.profiles,output.bowties,output.sams):
            shell("metaphlan2.py {f},{r} {p} --bowtie2out {b} --samout {s} --input_type fastq --nproc 4")


rule metaphlan2_merge:
    input:
        expand(os.path.join(output_folder, "metaphlan2/main/{samples}_profile.txt"), samples = SAMPLES)
    output:
        os.path.join(output_folder, "metaphlan2/merged/merged_abundance_table.txt")
    shell:
        "merge_metaphlan_tables.py {input} > {output}"


rule metaphlan2_report:
    input:
        hclust = os.path.join(output_folder, "metaphlan2/merged/abundance_heatmap_species.png"),
        graphlan = os.path.join(output_folder, "metaphlan2/merged/merged_abundance.png")
    output:
        "testing/metaphlan2/metaphlan2_report.html"
    run:
        from snakemake.utils import report
        report("""
        MetaPhlAn2 works!!!
        """, output[0])

# set options with "--config key1=value1 key2=value2"

import os

if "input_folder" in config:
    input_folder = config["input_folder"]
else:
    input_folder = "./"

if "output_folder" in config:
    output_folder = config["output_folder"]
else:
    output_folder = "./"


SAMPLES, = glob_wildcards(os.path.join(input_folder, "{samples}_R1_001.fastq.gz"))

rule all:
    input:
        os.path.join(input_folder, "report.html")


rule kneaddata_filter:
    input:
        fwd = expand(os.path.join(input_folder, "{samples}_R1_001.fastq.gz"), samples = SAMPLES),
        rev = expand(os.path.join(input_folder, "{samples}_R2_001.fastq.gz"), samples = SAMPLES),
        db = "testing/data/Homo_sapiens_Bowtie2_v0.1"
    output:
        folder = os.path.join(out_folder, "kneaddata/kneaddata_output/",
        fwd = expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_1.fastq"), samples = SAMPLES),
        rev = expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_2.fastq"), samples = SAMPLES)
    run:
        for f,r in zip(input.fwd,input.rev):
            shell("kneaddata --input {f} --input {r} --reference-db {input.db} --output {output.folder}")

        # def f():
        #   shell("ls testing/kneaddata/kneaddata_output/*kneaddata_paired* | while read F; do mv $F $( echo ${F} | sed 's/_kneaddata_paired/_paired/') ; done")
        #   mvfromF = '_'.join([wildcards.sample.replace('_',''), wildcards.read.replace('_kneaddata_paired',''), '_paired_1'])
        #   #shell('mv {mvfromF} {output.fwd}')
        #   #shell('mv {mvfromR} {output.rev')


rule kneaddata_counts:
    input:
        os.path.join(output_folder, "kneaddata/kneaddata_output/")
    output:
        os.path.join(output_folder, "kneaddata/kneaddata_read_counts.txt")
    shell:
        "kneaddata_read_count_table --input {input} --output {output}"


# rule kneaddata_clean:
#    input:
#        fwd = expand("testing/kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_1.fastq", samples = SAMPLES),
#        rev = expand("testing/kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_2.fastq", samples = SAMPLES),
#        # folder = "testing/kneaddata/kneaddata_output/"
#    output:
#        fwd = expand("testing/kneaddata/kneaddata_output/{samples}_R1_001_paired_1.fastq", samples = SAMPLES),
#        rev = expand("testing/kneaddata/kneaddata_output/{samples}_R1_001_paired_2.fastq", samples = SAMPLES)
#    run:
#        # shell("rm testing/kneaddata/kneaddata_output/*.fastq.bowtie2out.txt")
#        for f,r,x,y in zip(input.fwd,input.rev,output.fwd,output.rev):
#            shell("mv {f} {x}")
#            shell("mv {r} {y}")


rule kneaddata_report:
    input:
        filter = os.path.join(output_folder, "kneaddata/kneaddata_output/"),
        counts = os.path.join(output_folder, "kneaddata/kneaddata_read_counts.txt")
    output:
        "testing/kneaddata/kneaddata_report.html"
    run:
        from snakemake.utils import report
        report("""
        KneadData works!!!
        """, output[0])


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


rule hclust2_species:
    input:
        os.path.join(output_folder, "metaphlan2/merged/merged_abundance_table.txt")
    output:
        os.path.join(output_folder, "metaphlan2/merged/merged_abundance_table_species.txt")
    shell:
        "grep -E \"(s__)|(^ID)\" {input} | grep -v \"t__\" | sed 's/^.*s__//g' > {output}"


rule hclust2_heatmap:
    input:
        os.path.join(output_folder, "metaphlan2/merged/merged_abundance_table_species.txt")
    output:
        os.path.join(output_folder, "metaphlan2/merged/abundance_heatmap_species.png")
    shell:
        "hclust2.py -i {input} -o {output} --ftop 25 \
        --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 \
        --slabel_size 6 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300"


rule graphlan_input:
    input:
        os.path.join(output_folder, "metaphlan2/merged/merged_abundance_table.txt")
    output:
        tree = os.path.join(output_folder, "metaphlan2/merged/merged_abundance.tree.txt"),
        annot = os.path.join(output_folder, "metaphlan2/merged/merged_abundance.annot.txt")
    shell:
        "export2graphlan.py --skip_rows 1,2 \
        -i {input} --tree {output.tree} --annotation {output.annot} \
        --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 --annotations 5,6 \
        --external_annotations 7 --min_clade_size 1"


rule graphlan_cladogram_1:
    input:
        tree = os.path.join(output_folder, "metaphlan2/merged/merged_abundance.tree.txt"),
        annot = os.path.join(output_folder, "metaphlan2/merged/merged_abundance.annot.txt")
    output:
        os.path.join(output_folder, "metaphlan2/merged/merged_abundance.xml")
    shell:
        "graphlan_annotate.py --annot {input.annot} \
        {input.tree} {output}"


rule graphlan_cladogram_2:
    input:
        os.path.join(output_folder, "metaphlan2/merged/merged_abundance.xml")
    output:
        os.path.join(output_folder, "metaphlan2/merged/merged_abundance.png")
    shell:
        "graphlan.py --dpi 300 {input} \
        {output} --external_legends"


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


rule strainphlan_markers:
    input:
        expand("/Users/laurentso/Desktop/repos/echo/workflow/testing/metaphlan2/main/{samples}.sam.bz2", samples = SAMPLES)
    output:
        expand("testing/strainphlan/main/{samples}.markers", samples = SAMPLES)
    run:
        for i in input:
            shell("sample2markers.py --ifn_samples {i} --input_type sam --output_dir . \
            --samtools_exe /Users/laurentso/miniconda3/pkgs/samtools-0.1.19/bin/samtools \
            --bcftools_exe /Users/laurentso/miniconda3/pkgs/samtools-0.1.19/bin/bcftools")


# rule strainphlan_clades:
#    input:
#        expand("testing/strainphlan/main/{samples}.markers", samples = SAMPLES)
#    output:
#        "testing/strainphlan/merged/clades.txt"
#    run:
#        for i in input:
#            shell("strainphlan.py --ifn_samples {i} --output_dir . --print_clades_only > {output} --nprocs_main")


rule humann2_prep:
    input:
        fwd = expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_1.fastq"), samples = SAMPLES),
        rev = expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_2.fastq"), samples = SAMPLES)
    output:
        expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}.fastq"), samples = SAMPLES)
    run:
        for f,r,o in zip(input.fwd,input.rev,output):
            shell("cat {f} {r} > {o}")


rule humann2_output:
    input:
        expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}.fastq"), samples = SAMPLES)
    output:
        folder = "testing/humann2/main/",
        samples = expand("testing/humann2/main/{samples}_genefamilies.tsv", samples = SAMPLES),
        path = expand("testing/humann2/main/{samples}_pathabundance.tsv", samples = SAMPLES)
    run:
        for i in zip(input):
            shell("humann2 --input {i} --output {output.folder} --metaphlan ~/Documents/metaphlan2/metaphlan2 \
            --nucleotide-database ~/Desktop/repos/echo/workflow/testing/data/humann2_database_downloads/chocophlan \
            --protein-database ~/Desktop/repos/echo/workflow/testing/data/humann2_database_downloads/uniref \
            --diamond ~/Documents/diamond-linux64/")


rule humann2_rename:
    input:
        expand(os.path.join(output_folder, "humann2/main/{samples}_genefamilies.tsv"), samples = SAMPLES)
    output:
        expand(os.path.join(output_folder, "humann2/main/{samples}_genefamilies-names.tsv"), samples = SAMPLES)
    run:
        for i,o in zip(input,output):
            shell("humann2_rename_table --input {i} --output {o} --names uniref90")

rule humann2_regroup_1:
    input:
        samples = expand(os.path.join(output_folder, "humann2/main/{samples}_genefamilies.tsv"), samples = SAMPLES)
    output:
        expand(os.path.join(output_folder, "humann2/regroup/{samples}_ecs.tsv"), samples = SAMPLES)
    run:
        for i,o in zip(input.samples,output):
            shell("humann2_regroup_table --input {i} --output {o} --groups uniref90_rxn")


rule humann2_relab_1:
    input:
        # folder = "testing/humann2/main/",
        gf = expand(os.path.join(output_folder, "humann2/main/{samples}_genefamilies.tsv"), samples = SAMPLES),
        path = expand(os.path.join(output_folder, "humann2/main/{samples}_pathabundance.tsv"), samples = SAMPLES),
        ec = expand(os.path.join(output_folder, "humann2/regroup/{samples}_ecs.tsv"), samples = SAMPLES)
    output:
        gf = expand(os.path.join(output_folder, "humann2/relab/{samples}_genefamilies_relab.tsv"), samples = SAMPLES),
        path = expand(os.path.join(output_folder, "humann2/relab/{samples}_pathabundance_relab.tsv"), samples = SAMPLES),
        ec = expand(os.path.join(output_folder, "humann2/relab/{samples}_ecs_relab.tsv"), samples = SAMPLES)
    run:
        for g,p,e,x,y,z in zip(input.gf,input.path,input.ec,output.gf,output.path,output.ec):
            shell("humann2_renorm_table -i {g} -o {x} -u relab")
            shell("humann2_renorm_table -i {p} -o {y} -u relab")
            shell("humann2_renorm_table -i {e} -o {z} -u relab")


rule humann2_join:
    input:
        one = os.path.join(output_folder, "humann2/main"),
        two = os.path.join(output_folder, "humann2/regroup"),
        regroup = expand(os.path.join(output_folder, "humann2/regroup/{samples}_ecs.tsv"), samples = SAMPLES)
    output:
        gf = os.path.join(output_folder, "humann2/merged/genefamilies.tsv"),
        path = os.path.join(output_folder, "humann2/merged/pathabundance.tsv"),
        ec = os.path.join(output_folder, "humann2/merged/ecs.tsv")
    run:
        shell("humann2_join_tables -i {input.one} -o {output.gf} --file_name genefamilies")
        shell("humann2_join_tables -i {input.one} -o {output.path} --file_name pathabundance")
        shell("humann2_join_tables -i {input.two} -o {output.ec} --file_name ecs")


rule humann2_relab_2:
    input:
        gf = expand(os.path.join(output_folder, "humann2/relab/{samples}_genefamilies_relab.tsv"), samples = SAMPLES),
        path = expand(os.path.join(output_folder, "humann2/relab/{samples}_pathabundance_relab.tsv"), samples = SAMPLES),
        ec = expand(os.path.join(output_folder, "humann2/relab/{samples}_ecs_relab.tsv"), samples = SAMPLES)
    output:
        gf = os.path.join(output_folder, "humann2/merged/genefamilies_relab.tsv"),
        path = os.path.join(output_folder, "humann2/merged/pathabundance_relab.tsv"),
        ec = os.path.join(output_folder, "humann2/merged/ecs_relab.tsv")
    run:
        shell("humann2_join_tables -i testing/humann2/relab -o {output.gf} --file_name genefamilies_relab")
        shell("humann2_join_tables -i testing/humann2/relab -o {output.path} --file_name pathabundance_relab")
        shell("humann2_join_tables -i testing/humann2/relab -o {output.ec} --file_name ecs_relab")


rule humann2_report:
    input:
        gf = os.path.join(output_folder, "humann2/merged/genefamilies_relab.tsv"),
        path = os.path.join(output_folder, "humann2/merged/pathabundance_relab.tsv"),
        ec = os.path.join(output_folder, "humann2/merged/ecs_relab.tsv")
    output:
        os.path.join(output_folder, "humann2/humann2_report.html")
    run:
        from snakemake.utils import report
        report("""
        HUMAnN2 works!!!
        """, output[0])

rule report:
    input:
        kneaddata = os.path.join(output_folder, "kneaddata/kneaddata_report.html"),
        metaphlan2 = os.path.join(output_folder, "metaphlan2/metaphlan2_report.html"),
        humann2 = os.path.join(output_folder, "humann2/humann2_report.html")
    output:
        os.path.join(output_folder, "report.html")
    run:
        from snakemake.utils import report
        report("""
        Pipeline works!!!
        """, output[0])

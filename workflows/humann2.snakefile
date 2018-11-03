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
        samples = expand(os.path.join(output_folder, "humann2/main/{samples}_genefamilies.tsv"), samples = SAMPLES),
        path = expand(os.path.join(output_folder, "humann2/main/{samples}_pathabundance.tsv"), samples = SAMPLES)
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
        expand(os.path.join(output_folder, "humann2/main/{samples}_genefamilies.tsv"), samples = SAMPLES)
    output:
        expand(os.path.join(output_folder, "humann2/regroup/{samples}_ecs.tsv"), samples = SAMPLES)
    run:
        for i,o in zip(input,output):
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

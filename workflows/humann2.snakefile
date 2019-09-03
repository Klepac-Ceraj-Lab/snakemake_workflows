##############
# Per Sample #
##############

rule humann2:
    input:
        seq = os.path.join(kneadfolder, "{sample}_merged.fastq"),
        tax = os.path.join(metaphlanfolder, "main", "{sample}_profile.tsv")
    output:
        samples = os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv"),
        path = os.path.join(humannfolder, "main", "{sample}_pathabundance.tsv")
    run:
        # TODO: get threads from settings
        shell("humann2 --input {{input.seq}} --taxonomic-profile {{input.tax}} --output {} --threads 8 --remove-temp-output --search-mode uniref90 --output-basename {{wildcards.sample}}".format(
            os.path.join(humannfolder, "main")))


rule humann2_regroup_ecs:
    input: os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv")
    output: os.path.join(humannfolder, "regroup", "{sample}_ecs.tsv")
    run: shell("humann2_regroup_table --input {input} --output {output} --groups uniref90_level4ec")

rule humann2_regroup_kos:
    input: os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv")
    output: os.path.join(humannfolder, "regroup", "{sample}_kos.tsv")
    run: shell("humann2_regroup_table --input {input} --output {output} --groups uniref90_ko")

rule humann2_renorm_gf:
    input: os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_genefamilies_relab.tsv")
    run: shell("humann2_renorm_table -i {input} -o {output} -u relab")

rule humann2_renorm_ecs:
    input: os.path.join(humannfolder, "regroup", "{sample}_ecs.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_ecs_relab.tsv")
    run: shell("humann2_renorm_table -i {input} -o {output} -u relab")

rule humann2_renorm_kos:
    input: os.path.join(humannfolder, "regroup", "{sample}_kos.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_kos_relab.tsv")
    run: shell("humann2_renorm_table -i {input} -o {output} -u relab")

rule humann2_renorm_paths:
    input: os.path.join(humannfolder, "main", "{sample}_pathabundance.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_pathabundance_relab.tsv")
    run: shell("humann2_renorm_table -i {input} -o {output} -u relab")


rule humann2_rename_gf:
    input: os.path.join(humannfolder, "relab", "{samples}_genefamilies_relab.tsv")
    output: os.path.join(humannfolder, "main", "{samples}_genefamilies_names.tsv")
    run:
        shell("humann2_rename_table --input {input} --output {output} --names uniref90")

rule humann2_rename_ecs:
    input: os.path.join(humannfolder, "relab", "{samples}_ecs_relab.tsv")
    output: os.path.join(humannfolder, "main", "{samples}_ecs_names.tsv")
    run:
        shell("humann2_rename_table --input {input} --output {output} --names uniref90")

rule humann2_rename_kos:
    input: os.path.join(humannfolder, "relab", "{samples}_kos_relab.tsv")
    output: os.path.join(humannfolder, "main", "{samples}_kos_names.tsv")
    run:
        shell("humann2_rename_table --input {input} --output {output} --names uniref90")


rule humann2_report:
    input:
        gf = expand(os.path.join(humannfolder, "main", "{sample}_genefamilies_names.tsv"), sample=samples),
        path = expand(os.path.join(humannfolder, "relab", "{sample}_pathabundance_relab.tsv"), sample=samples),
        ec = expand(os.path.join(humannfolder, "main", "{sample}_ecs_names.tsv"), sample=samples),
        ko = expand(os.path.join(humannfolder, "main", "{sample}_kos_names.tsv"), sample=samples)
    output:
        os.path.join(humannfolder, "humann2_report.html")
    run:
        from snakemake.utils import report
        report("""
        HUMAnN2 works!!!
        """, output[0])

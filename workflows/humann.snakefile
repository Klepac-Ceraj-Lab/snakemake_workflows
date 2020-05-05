##############
# Per Sample #
##############

rule humann:
    input:
        seq = os.path.join(kneadfolder, "{sample}_merged.fastq"),
        tax = os.path.join(metaphlanfolder, "{sample}_profile.tsv")
    output:
        samples = os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv"),
        path = os.path.join(humannfolder, "main", "{sample}_pathabundance.tsv")
    run:
        # TODO: get threads from settings
        shell("humann --input {{input.seq}} --taxonomic-profile {{input.tax}} --output {} --threads 8 --remove-temp-output --search-mode uniref90 --output-basename {{wildcards.sample}} --metaphlan-options '-x mpa_v296_CHOCOPhlAn_201901'".format(
            os.path.join(humannfolder, "main")))


rule humann_regroup_ecs:
    input: os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv")
    output: os.path.join(humannfolder, "main", "{sample}_ecs.tsv")
    run: shell("humann_regroup_table --input {input} --output {output} --groups uniref90_level4ec")

rule humann_regroup_kos:
    input: os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv")
    output: os.path.join(humannfolder, "main", "{sample}_kos.tsv")
    run: shell("humann_regroup_table --input {input} --output {output} --groups uniref90_ko")

rule humann_regroup_pfams:
    input: os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv")
    output: os.path.join(humannfolder, "main", "{sample}_pfams.tsv")
    run: shell("humann_regroup_table --input {input} --output {output} --groups uniref90_pfam")

rule humann_renorm_gf:
    input: os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_genefamilies_relab.tsv")
    run: shell("humann_renorm_table -i {input} -o {output} -u relab")

rule humann_renorm_ecs:
    input: os.path.join(humannfolder, "main", "{sample}_ecs.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_ecs_relab.tsv")
    run: shell("humann_renorm_table -i {input} -o {output} -u relab")

rule humann_renorm_kos:
    input: os.path.join(humannfolder, "main", "{sample}_kos.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_kos_relab.tsv")
    run: shell("humann_renorm_table -i {input} -o {output} -u relab")

rule humann_renorm_pfams:
    input: os.path.join(humannfolder, "main", "{sample}_pfams.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_pfams_relab.tsv")
    run: shell("humann_renorm_table -i {input} -o {output} -u relab")

rule humann_renorm_paths:
    input: os.path.join(humannfolder, "main", "{sample}_pathabundance.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_pathabundance_relab.tsv")
    run: shell("humann_renorm_table -i {input} -o {output} -u relab")


rule humann_rename_ecs:
    input: os.path.join(humannfolder, "relab", "{samples}_ecs_relab.tsv")
    output: os.path.join(humannfolder, "names", "{samples}_ecs_names_relab.tsv")
    run: shell("humann_rename_table --input {input} --output {output} --names ec")

rule humann_rename_kos:
    input: os.path.join(humannfolder, "relab", "{samples}_kos_relab.tsv")
    output: os.path.join(humannfolder, "names", "{samples}_kos_names_relab.tsv")
    run: shell("humann_rename_table --input {input} --output {output} --names kegg-orthology")

rule humann_rename_pfams:
    input: os.path.join(humannfolder, "relab", "{samples}_pfams_relab.tsv")
    output: os.path.join(humannfolder, "names", "{samples}_pfams_names_relab.tsv")
    run: shell("humann_rename_table --input {input} --output {output} --names pfam")

rule humann_report:
    input:
        gf = expand(os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv"), sample=samples),
        path = expand(os.path.join(humannfolder, "relab", "{sample}_pathabundance_relab.tsv"), sample=samples),
        ko = expand(os.path.join(humannfolder, "names", "{sample}_kos_names_relab.tsv"), sample=samples),
        ec = expand(os.path.join(humannfolder, "names", "{sample}_ecs_names_relab.tsv"), sample=samples),
        pfam = expand(os.path.join(humannfolder, "names", "{sample}_pfams_names_relab.tsv"), sample=samples)
    output:
        os.path.join(humannfolder, "humann_report.html")
    run:
        from snakemake.utils import report
        report("""
        humann works!!!
        """, output[0])

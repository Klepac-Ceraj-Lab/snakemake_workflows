humannfolder = os.path.join(output_folder, "humann")

##############
# Per Sample #
##############

rule humann:
    input:
        seq = os.path.join(kneadfolder, "{sample}_kneaddata.fastq"),
        tax = os.path.join(metaphlanfolder, "{sample}_profile.tsv")
    output:
        samples = os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv"),
        path = os.path.join(humannfolder, "main", "{sample}_pathabundance.tsv")
    run:
        shell("humann --input {{input.seq}} --taxonomic-profile {{input.tax}} --output {} --threads {cluster['processors']} --remove-temp-output --search-mode uniref90 --output-basename {{wildcards.sample}}".format(
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

rule humann_rename_ecs:
    input: os.path.join(humannfolder, "main", "{samples}_ecs.tsv")
    output: os.path.join(humannfolder, "names", "{samples}_ecs_names.tsv")
    run: shell("humann_rename_table --input {input} --output {output} --names ec")

rule humann_rename_kos:
    input: os.path.join(humannfolder, "main", "{samples}_kos.tsv")
    output: os.path.join(humannfolder, "names", "{samples}_kos_names.tsv")
    run: shell("humann_rename_table --input {input} --output {output} --names kegg-orthology")

rule humann_rename_pfams:
    input: os.path.join(humannfolder, "main", "{samples}_pfams.tsv")
    output: os.path.join(humannfolder, "names", "{samples}_pfams_names.tsv")
    run: shell("humann_rename_table --input {input} --output {output} --names pfam")

rule humann_report:
    input:
        gf = expand(os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv"), sample=samples),
        ko = expand(os.path.join(humannfolder, "names", "{sample}_kos_names.tsv"), sample=samples),
        ec = expand(os.path.join(humannfolder, "names", "{sample}_ecs_names.tsv"), sample=samples),
        pfam = expand(os.path.join(humannfolder, "names", "{sample}_pfams_names.tsv"), sample=samples)
    output:
        os.path.join(humannfolder, "humann_report.html")
    run:
        from snakemake.utils import report
        report("""
        humann works!!!
        """, output[0])

##############
# Per Sample #
##############

rule humann2:
    input:
        os.path.join(metaphlanfolder, "main", "{sample}.sam")
    output:
        samples = os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv"),
        path = os.path.join(humannfolder, "main", "{sample}_pathabundance.tsv")
    run:
        # TODO: get threads from settings
        shell("humann2 --input {{input}} --output {} --threads 8 --nucleotide-database {} --protein-database {} --input-format sam".format(
            os.path.join(humannfolder, "main"),
            config["databases"]["chocophlan"],
            config["databases"]["uniref"]))


rule humann2_regroup_ecs:
    input: os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv")
    output: os.path.join(humannfolder, "regroup", "{sample}_ecs.tsv")
    run: shell("humann2_regroup_table --input {input} --output {output} --groups uniref90_rxn")

rule humann2_renorm_gf:
    input: os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_genefamilies_relab.tsv")
    run: shell("humann2_renorm_table -i {input} -o {output} -u relab")

rule humann2_renorm_ecs:
    input: os.path.join(humannfolder, "regroup", "{sample}_ecs.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_ecs_relab.tsv")
    run: shell("humann2_renorm_table -i {input} -o {output} -u relab")

rule humann2_renorm_paths:
    input: os.path.join(humannfolder, "main", "{sample}_pathabundance.tsv")
    output: os.path.join(humannfolder, "relab", "{sample}_pathabundance_relab.tsv")
    run: shell("humann2_renorm_table -i {input} -o {output} -u relab")


# rule humann2_rename_gf:
#     input: os.path.join(humannfolder, "main", "{samples}_genefamilies.tsv")
#     output: os.path.join(humannfolder, "main", "{samples}_genefamilies_names.tsv"),
#     run:
#         shell("humann2_rename_table --input {input} --output {output} --names uniref90")



###############
# All samples #
###############

rule humann2_merge_gf:
    input: expand(os.path.join(humannfolder, "main", "{sample}_genefamilies.tsv"), sample = samples)
    output: os.path.join(humannfolder, "merged", "genefamilies.tsv")
    run: shell("humann2_join_tables -i {input} -o {output} --file_name genefamilies")

rule humann2_merge_ecs:
    input: expand(os.path.join(humannfolder, "regroup", "{sample}_ecs.tsv"), sample = samples)
    output: os.path.join(humannfolder, "merged", "ecs.tsv")
    run: shell("humann2_join_tables -i {input} -o {output} --file_name ecs")

rule humann2_merge_paths:
    input: expand(os.path.join(humannfolder, "main", "{sample}_pathabundance.tsv"), sample = samples)
    output: os.path.join(humannfolder, "merged", "pathabundance.tsv")
    run: shell("humann2_join_tables -i {input} -o {output} --file_name pathabundance")

rule humann2_merge_gf_relab:
    input: expand(os.path.join(humannfolder, "relab", "{sample}_genefamilies_relab.tsv"), sample = samples)
    output: os.path.join(humannfolder, "merged", "genefamilies_relab.tsv")
    run: shell("humann2_join_tables -i {} -o {{output}} --file_name genefamilies_relab".format(os.path.join(humannfolder, "merged")))

rule humann2_merge_ecs_relab:
    input: expand(os.path.join(humannfolder, "relab", "{sample}_ecs_relab.tsv"), sample = samples)
    output: os.path.join(humannfolder, "merged", "ecs_relab.tsv")
    run: shell("humann2_join_tables -i {} -o {{output}} --file_name ecs_relab".format(os.path.join(humannfolder, "merged")))

rule humann2_merge_paths_relab:
    input: expand(os.path.join(humannfolder, "relab", "{sample}_pathabundance_relab.tsv"), sample = samples)
    output: os.path.join(humannfolder, "merged", "pathabundance_relab.tsv")
    run: shell("humann2_join_tables -i {} -o {{output}} --file_name pathabundance_relab".format(os.path.join(humannfolder, "merged")))


rule humann2_report:
    input:
        gf = os.path.join(humannfolder, "merged", "genefamilies_relab.tsv"),
        path = os.path.join(humannfolder, "merged", "pathabundance_relab.tsv"),
        ec = os.path.join(humannfolder, "merged", "ecs_relab.tsv")
    output:
        os.path.join(humannfolder, "humann2_report.html")
    run:
        from snakemake.utils import report
        report("""
        HUMAnN2 works!!!
        """, output[0])

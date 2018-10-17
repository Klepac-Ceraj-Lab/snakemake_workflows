## Not functional

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

rule hclust_graphlan_report:
    input:
        hclust = os.path.join(output_folder, "metaphlan2/merged/abundance_heatmap_species.png"),
        graphlan = os.path.join(output_folder, "metaphlan2/merged/merged_abundance.png")
    output:
        "testing/metaphlan2/hclust_graphlan_report.html"
    run:
        from snakemake.utils import report
        report("""
        Hclust and GraphPhlAn work!!!
        """, output[0])

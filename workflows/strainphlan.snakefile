## Not functional


rule strainphlan_markers:
    input:
        expand("/Users/laurentso/Desktop/repos/echo/workflow/testing/metaphlan/main/{samples}.sam.bz2", samples = SAMPLES)
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

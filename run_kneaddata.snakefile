kneadfolder = os.path.join(output_folder, "kneaddata")
include: "workflows/kneaddata.snakefile"

rule all:
    input:
        os.path.join(output_folder, "report.html")


rule report:
    input:
        kneaddata = os.path.join(kneadfolder, "kneaddata_report.html"),
    output:
        os.path.join(output_folder, "report.html")
    run:
        from snakemake.utils import report
        report("""
        Pipeline works!!!
        """, output[0])

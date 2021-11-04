include: "workflows/kneaddata.snakefile"
include: "workflows/metaphlan.snakefile"
include: "workflows/humann.snakefile"
include: "workflows/cleanup.snakefile"

rule report:
    input:
        cleanup = os.path.join(output_folder, "cleanup_report.html"),
    output:
        os.path.join(output_folder, "report.html")
    run:
        from snakemake.utils import report
        report("""
        Pipeline works!!!
        """, output[0])
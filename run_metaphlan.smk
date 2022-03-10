include: "workflows/kneaddata.smk"
include: "workflows/metaphlan.smk"

rule report:
    input:
        kneaddata = os.path.join(kneadfolder, "kneaddata_report.html"),
        metaphlan = os.path.join(metaphlanfolder, "metaphlan_report.html"),
    output:
        os.path.join(output_folder, "report.html")
    run:
        from snakemake.utils import report
        report("""
        Pipeline works!!!
        """, output[0])

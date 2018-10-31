rule kneaddata_filter:
    input:
        fwd = expand(os.path.join(input_folder, "{samples}_R1_001.fastq.gz"), samples = SAMPLES),
        rev = expand(os.path.join(input_folder, "{samples}_R2_001.fastq.gz"), samples = SAMPLES),
        db = config["databases"]["human_sequences"]
    output:
        folder = os.path.join(out_folder, "kneaddata/kneaddata_output/"),
        fwd = expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_1.fastq"), samples = SAMPLES),
        rev = expand(os.path.join(output_folder, "kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_2.fastq"), samples = SAMPLES)
    run:
        for f,r in zip(input.fwd,input.rev):
            shell("kneaddata --input {f} --input {r} --reference-db {input.db} --output {output.folder}")

        # def f():
        #   shell("ls testing/kneaddata/kneaddata_output/*kneaddata_paired* | while read F; do mv $F $( echo ${F} | sed 's/_kneaddata_paired/_paired/') ; done")
        #   mvfromF = '_'.join([wildcards.sample.replace('_',''), wildcards.read.replace('_kneaddata_paired',''), '_paired_1'])
        #   #shell('mv {mvfromF} {output.fwd}')
        #   #shell('mv {mvfromR} {output.rev')


rule kneaddata_counts:
    input:
        os.path.join(output_folder, "kneaddata/kneaddata_output/")
    output:
        os.path.join(output_folder, "kneaddata/kneaddata_read_counts.txt")
    shell:
        "kneaddata_read_count_table --input {input} --output {output}"


# rule kneaddata_clean:
#    input:
#        fwd = expand("testing/kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_1.fastq", samples = SAMPLES),
#        rev = expand("testing/kneaddata/kneaddata_output/{samples}_R1_001_kneaddata_paired_2.fastq", samples = SAMPLES),
#        # folder = "testing/kneaddata/kneaddata_output/"
#    output:
#        fwd = expand("testing/kneaddata/kneaddata_output/{samples}_R1_001_paired_1.fastq", samples = SAMPLES),
#        rev = expand("testing/kneaddata/kneaddata_output/{samples}_R1_001_paired_2.fastq", samples = SAMPLES)
#    run:
#        # shell("rm testing/kneaddata/kneaddata_output/*.fastq.bowtie2out.txt")
#        for f,r,x,y in zip(input.fwd,input.rev,output.fwd,output.rev):
#            shell("mv {f} {x}")
#            shell("mv {r} {y}")


rule kneaddata_report:
    input:
        filter = os.path.join(output_folder, "kneaddata/kneaddata_output/"),
        counts = os.path.join(output_folder, "kneaddata/kneaddata_read_counts.txt")
    output:
        "testing/kneaddata/kneaddata_report.html"
    run:
        from snakemake.utils import report
        report("""
        KneadData works!!!
        """, output[0])

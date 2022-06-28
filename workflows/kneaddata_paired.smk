rule kneaddata_cat_pair1:
    input: expand(os.path.join(input_folder, "{{sample}}_{lane}_R1_001.fastq.gz"), lane=lanes)
    output: os.path.join(input_folder, "{sample}_1.fastq.gz")
    run:
        shell("cat {input} > {output}")
rule kneaddata_cat_pair2:
    input: expand(os.path.join(input_folder, "{{sample}}_{lane}_R2_001.fastq.gz"), lane=lanes)
    output: os.path.join(input_folder, "{sample}_2.fastq.gz")
    run:
        shell("cat {input} > {output}")

rule kneaddata:
    input:
        fwd = os.path.join(input_folder, "{sample}_1.fastq.gz"),
        rev = os.path.join(input_folder, "{sample}_2.fastq.gz"),
    output:
        fwd = temp(os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq")),
        rev = temp(os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq")),
        log = os.path.join(kneadfolder, "{sample}_kneaddata.log")
    run:
        shell("kneaddata --input {{input.fwd}} --input {{input.rev}} --reference-db /hg37 --output {} --output-prefix {{wildcards.sample}}_kneaddata --trimmomatic /opt/conda/share/trimmomatic".format(kneadfolder))
    
rule compressdata1:
    input: 
        fwd= os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"),
    output:
        rev= os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq.gz")
    run: 
        shell("gzip -c {input} > {output}")

rule compressdata2:
    input: 
        fwd= os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq"),
    output:
        rev= os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq.gz")
    run: 
        shell("gzip -c {input} > {output}")

####trying for loop
#checkpoint further_compress:
#        directory = os.fsencode(kneadfolder)
#        IDs = []
#    for file in os.listdir(directory):
#        filename = os.fsdecode(file)
#        if filename.endswith(".fastq"):
#            IDs.append(filename)

#    for id in IDs:
#	    input: '{id}'
#	    output: '{id}.gz'
#	    run: shell("gzip -c {input} > {output}")

 
####trying with glob
#IDs, = glob_wildcards("kneadfolder/{id}.fastq")
#print("these are the IDS", IDs)

checkpoint all_:
   output:
    expand("kneadfolder/{sample}.fastq.gz",
    sample = samples)

def aggregate_input(wildcards):
    samples = glob_wildcards("kneadfolder/{sample}.fastq")
    checkpoints.all_.get(sample=wildcards.sample)
    return ("kneadfolder/{sample}.fastq")

rule further_compress:
   input: aggregate_input
   output: ("kneadfolder/{sample}.fastq.gz")
   run: shell("cp {input} {output}")


#rule further_compress:
	#input: '{sample}_{suffix}.fastq'
	#output: '{sample}_{suffix}.fastq.gz'
	#shell: ("gzip {input}")

rule metaphlan_cat:
    input:
        fwd = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"),
        rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq"),
    output: temp(os.path.join(kneadfolder, "{sample}_kneaddata.fastq"))
    run:
        shell("cat {input} > {output}")

rule kneaddata_counts:
    input: expand(os.path.join(kneadfolder, "{sample}_kneaddata.log"), sample = samples)
    output: os.path.join(kneadfolder, "kneaddata_read_counts.tsv")
    shell:
        "kneaddata_read_count_table --input {} --output {{output}}".format(kneadfolder)

rule kneaddata_report:
    input:
        counts = os.path.join(kneadfolder, "kneaddata_read_counts.tsv"),
        fwd = expand(os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq.gz"), sample = samples),
        rev = expand(os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq.gz"), sample = samples)
    output:
        os.path.join(kneadfolder, "kneaddata_report.html")
    run:
        from snakemake.utils import report
        report("""
        KneadData works!!!
        """, output[0])

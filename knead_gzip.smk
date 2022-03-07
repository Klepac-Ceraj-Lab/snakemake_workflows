import os, glob

configfile: "config.yaml"

include: "setup/directories.snakefile"
include: config["samples"]

kneadfolder = os.path.join(output_folder, "kneaddata")

rule all:
    input:
        os.path.join(output_folder, "knead_zipped.html")

rule zipped_report:
    input:
        paired1 = expand(os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq.gz"), sample=samples),
        repeats1 = expand(os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.1.fastq.gz"), sample=samples),
        repeats_unmatched1 = expand(os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.unmatched.1.fastq.gz"), sample=samples),
        trimmed1 = expand(os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.1.fastq.gz"), sample=samples),
        trimmed_single1 = expand(os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.single.1.fastq.gz"), sample=samples),
        unmatched1 = expand(os.path.join(kneadfolder, "{sample}_kneaddata_unmatched_1.fastq.gz"), sample=samples),
        hg371 = expand(os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_paired_contam_1.fastq.gz"), sample=samples),
        hg37_unmatched1 = expand(os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_unmatched_1_contam.fastq.gz"), sample=samples),
        paired2 = expand(os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq.gz"), sample=samples),
        repeats2 = expand(os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.2.fastq.gz"), sample=samples),
        repeats_unmatched2 = expand(os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.unmatched.2.fastq.gz"), sample=samples),
        trimmed2 = expand(os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.2.fastq.gz"), sample=samples),
        trimmed_single2 = expand(os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.single.2.fastq.gz"), sample=samples),
        unmatched2 = expand(os.path.join(kneadfolder, "{sample}_kneaddata_unmatched_2.fastq.gz"), sample=samples),
        hg372 = expand(os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_paired_contam_2.fastq.gz"), sample=samples),
        hg37_unmatched2 = expand(os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_unmatched_2_contam.fastq.gz"), sample=samples),
    output:
        os.path.join(output_folder, "knead_zipped.html")
    run:
        from snakemake.utils import report
        report("""
        KneadData works!!!
        """, output[0])

rule kneaddata_pair1:
    input:
        paired = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"),
        repeats = os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.1.fastq"),
        repeats_unmatched = os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.unmatched.1.fastq"),
        trimmed = os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.1.fastq"),
        trimmed_single = os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.single.1.fastq"),
        unmatched = os.path.join(kneadfolder, "{sample}_kneaddata_unmatched_1.fastq"),
        hg37 = os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_paired_contam_1.fastq"),
        hg37_unmatched = os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_unmatched_1_contam.fastq")
    output:
        paired = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq.gz"),
        repeats = os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.1.fastq.gz"),
        repeats_unmatched = os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.unmatched.1.fastq.gz"),
        trimmed = os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.1.fastq.gz"),
        trimmed_single = os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.single.1.fastq.gz"),
        unmatched = os.path.join(kneadfolder, "{sample}_kneaddata_unmatched_1.fastq.gz"),
        hg37 = os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_paired_contam_1.fastq.gz"),
        hg37_unmatched = os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_unmatched_1_contam.fastq.gz")
    run:
        for v in input:
            shell("gzip -f {v} --rsyncable")

rule kneaddata_pair2:
    input:
        paired = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq"),
        repeats = os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.2.fastq"),
        repeats_unmatched = os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.unmatched.2.fastq"),
        trimmed = os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.2.fastq"),
        trimmed_single = os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.single.2.fastq"),
        unmatched = os.path.join(kneadfolder, "{sample}_kneaddata_unmatched_2.fastq"),
        hg37 = os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_paired_contam_2.fastq"),
        hg37_unmatched = os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_unmatched_2_contam.fastq")
    output:
        paired = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq.gz"),
        repeats = os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.2.fastq.gz"),
        repeats_unmatched = os.path.join(kneadfolder, "{sample}_kneaddata.repeats.removed.unmatched.2.fastq.gz"),
        trimmed = os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.2.fastq.gz"),
        trimmed_single = os.path.join(kneadfolder, "{sample}_kneaddata.trimmed.single.2.fastq.gz"),
        unmatched = os.path.join(kneadfolder, "{sample}_kneaddata_unmatched_2.fastq.gz"),
        hg37 = os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_paired_contam_2.fastq.gz"),
        hg37_unmatched = os.path.join(kneadfolder, "{sample}_kneaddata_hg37dec_v0.1_bowtie2_unmatched_2_contam.fastq.gz")
    run:
        for v in input:
            shell("gzip -f {v} --rsyncable")


